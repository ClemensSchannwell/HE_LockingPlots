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
    if args.region == "Kenzie":
        base_sample = get_he_events("HE02", args.region, threshold=1.0e12)
        grad_thresh = -2e11
        flux_thresh = 1.5e12
        figno = "04"
    else:
        base_sample = get_he_events("HE02", args.region)
        grad_thresh = -3e11
        flux_thresh = 2.25e12
        figno = "02"
    t_pert = [14, 10, 6, 2]
    smb_pert = [150, 100, 75, 50, 37.5, 9.375]
    runs = [["HE68", "HE69", "HE70", "HE71", "HE72", "HE73"],
            ["HE46", "HE61", "HE55", "HE64", "HE49", "HE52"],
            ["HE47", "HE62", "HE56", "HE65", "HE50", "HE53"],
            ["HE48", "HE63", "HE57", "HE66", "HE51", "HE54"]]
    ticks = np.array([0, 300, 600, 900, 1200])
    xticks = 2 * np.pi / 1500.0 * ticks
    labels = ['0', '300', '600', '900', '1200']
    m, n = 6, 4
    fig, ax = plt.subplots(m, n, subplot_kw={'projection': 'polar'},
                           gridspec_kw={'wspace': 0.15, 'hspace': 0.175,
                                        'top': 0.99, 'bottom': 0.01,
                                        'left': 0.03, 'right': 0.97},
                           figsize=(10.5, 14))
    for i_smb in range(len(runs[0])):
        for i_temp in range(len(runs)):
            fname = f'Data/IceVolume{args.region}_{runs[i_temp][i_smb]}.nc'
            print(fname)
            time, vol = read_data(fname)
            he_sample = get_he_events(runs[i_temp][i_smb], args.region, flux_thresh)
            sig_level = kuipers_test(he_sample, base_sample)
            plot_siglevel(sig_level, ax[i_smb][i_temp])
            surge_events = extract_surges(vol, thresh=grad_thresh)
            theta = 2 * np.pi / 1500.0 * time
            plot_he_lines(ax[i_smb][i_temp], theta, vol, surge_events)
            plot_do_cycle_polar(ax[i_smb][i_temp])
            set_axis_properties(args.region, ax[i_smb][i_temp],
                                xticks, labels)
            plot_panel_titles(ax[i_smb][i_temp], i_smb, i_temp, t_pert,
                              smb_pert)
            create_subplot_labels(ax[i_smb][i_temp], i_smb, i_temp)
    ax[0, 0].legend(loc=1, bbox_to_anchor=(1.35, 1.175))
    plt.savefig(f'ExFig{figno}_Polar_{args.region}.png', dpi=300,
                bbox_inches='tight')


if __name__ == '__main__':
    main()
