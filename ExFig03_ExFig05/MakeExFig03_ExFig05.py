#! /usr/bin/env python 

import argparse
import sys
sys.path.append("../util")
from utility_fcts import *


def create_subplot_labels(ax, row, col):
    label_mat = [['(a)', '(b)', '(c)', '(d)'], ['(e)', '(f)', '(g)', '(h)'],
                 ['(i)', '(j)', '(k)', '(l)'], ['(m)', '(n)', '(o)', '(p)'],
                 ['(q)', '(r)', '(s)', '(t)'], ['(u)', '(v)', '(w)', '(x)'],
                 ['(y)', '(z)', '(aa)', '(bb)']]
    ax.text(0.05, 0.95, label_mat[row][col], color='black',
            horizontalalignment='center', verticalalignment='center',
            transform=ax.transAxes, fontsize=14)


def plot_panel_titles(ax, i_smb, i_temp, t_pert, smb_pert):
    if i_smb == 0:
        ax.set_title("$\Delta$Temp = " + str(t_pert[i_temp]) + "°C",
                     fontsize=14, fontweight='bold')
    if i_temp == 0:
        smb_str = "$\Delta$SMB = " + str(round(smb_pert[i_smb] / 1000, 3)) + " m/yr"
        ax.text(-0.2, 0.5, smb_str, color="black", horizontalalignment='center',
                verticalalignment='center', transform=ax.transAxes, rotation=90,
                fontsize=14, fontweight='bold')


def set_axis_properties(region, ax, xticks, labels):
    if region == "Kenzie":
        ax.set_ylim(3.85e15, 4.95e15)
    else:
        ax.set_ylim(4e15, 5.9e15)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_yticklabels([])
    ax.yaxis.grid(False)
    ax.set_xticks(xticks)
    ax.set_xticklabels(labels)
    ax.xaxis.grid(True)
    ax.tick_params(pad=2)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-reg", "--region", choices=["Hudson", "Kenzie"],
                        help="Specify which plot you want to create",
                        required=True)
    args = parser.parse_args()
    clen = [1000, 2000]
    t_pert = [14, 10, 6, 2]
    smb_pert = [100, 50]
    for i_clen in clen:
        if args.region == "Kenzie":
            base_sample = get_he_events("HE02", args.region, threshold=1.0e12,
                                        clen=i_clen)
            grad_thresh = -2e11
            flux_thresh = 1.5e12
            figno = "05"
        else:
            base_sample = get_he_events("HE02", args.region, clen=i_clen)
            grad_thresh = -3e11
            flux_thresh = 2.25e12
            figno = "03"
        if i_clen == 1000:
            runs = [["HE81", "HE86"], ["HE78", "HE87"], ["HE79", "HE88"],
                    ["HE80", "HE89"]]
            ticks = np.array([0, 250, 500, 750])
            labels = ['0', '250', '500', '750']
        elif i_clen == 2000:
            runs = [["HE77", "HE82"], ["HE74", "HE83"], ["HE75", "HE84"],
                    ["HE76", "HE85"]]
            ticks = np.array([0, 500, 1000, 1500])
            labels = ['0', '500', '1000', '1500']

        xticks = 2 * np.pi / i_clen * ticks
        m, n = 4, 2
        fig, ax = plt.subplots(n, m, subplot_kw={'projection': 'polar'},
                               gridspec_kw={'wspace': 0.43, 'hspace': 0.175,
                                            'top': 0.99, 'bottom': 0.01,
                                            'left': 0.03, 'right': 0.97},
                               figsize=(12, 5))
        for i_smb in range(len(runs[0])):
            for i_temp in range(len(runs)):
                fname = f'Data/IceVolume{args.region}_{runs[i_temp][i_smb]}.nc'
                time, vol = read_data(fname)
                he_sample = get_he_events(runs[i_temp][i_smb], args.region,
                                          flux_thresh, clen=i_clen)
                sig_level = kuipers_test(he_sample, base_sample)
                plot_siglevel(sig_level, ax[i_smb][i_temp])
                surge_events = extract_surges(vol, thresh=grad_thresh)
                theta = 2 * np.pi / i_clen * time
                plot_he_lines(ax[i_smb][i_temp], theta, vol, surge_events)
                plot_do_cycle_polar(ax[i_smb][i_temp], clen=i_clen)
                set_axis_properties(args.region, ax[i_smb][i_temp],
                                    xticks, labels)
                plot_panel_titles(ax[i_smb][i_temp], i_smb, i_temp, t_pert,
                                  smb_pert)
                create_subplot_labels(ax[i_smb][i_temp], i_smb, i_temp)
        ax[0, 0].legend(loc=1, bbox_to_anchor=(1.25, 1.1))
        plt.savefig(f'ExFig{figno}_Polar_{args.region}_{i_clen}.png',
                    dpi=300, bbox_inches='tight')
    plt.show()


if __name__ == '__main__':
    main()
