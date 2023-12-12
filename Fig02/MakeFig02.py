#! /usr/bin/env python

import argparse
import sys
sys.path.append("../util")
from utility_fcts import *


def create_subplot_labels(ax, row, col):
    label_mat = [['(a)', '(b)'], ['(c)', '(d)'], ['(e)', '(f)']]
    ax.text(0.05, 0.95, label_mat[row][col], color='black',
            horizontalalignment='center', verticalalignment='center',
            transform=ax.transAxes, fontsize=17)


def main():
    # global base_path, region, t_pert, smb_pert
    parser = argparse.ArgumentParser()
    parser.add_argument("-reg", "--region", choices=["Hudson", "Kenzie"],
                        help="Specify which plot you want to create",
                        required=True)
    args = parser.parse_args()
    if args.region == "Kenzie":
        base_sample = get_he_events("HE02", args.region, threshold=1.0e12)
        grad_thresh = -2e11
        flux_thresh = 1.5e12
    else:
        base_sample = get_he_events("HE02", args.region)
        grad_thresh = -3e11
        flux_thresh = 2.25e12
    t_pert = [14, 10, 6, 2]
    smb_pert = [150, 100, 75, 50, 37.5, 9.375]
    runs = [["HE68", "HE69", "HE70", "HE71", "HE72", "HE73"],
            ["HE46", "HE61", "HE55", "HE64", "HE49", "HE52"],
            ["HE47", "HE62", "HE56", "HE65", "HE50", "HE53"],
            ["HE48", "HE63", "HE57", "HE66", "HE51", "HE54"]]
    ticks = np.array([0, 300, 600, 900, 1200])
    xticks = 2 * np.pi / 1500.0 * ticks
    labels = ['0', '300', '600', '900', '1200']
    m, n = 3, 2
    fig, ax = plt.subplots(m, n, subplot_kw={'projection': 'polar'},
                           gridspec_kw={'wspace': 0.15, 'hspace': 0.175,
                                        'top': 0.99, 'bottom': 0.01,
                                        'left': 0.03, 'right': 0.97},
                           figsize=(10.5, 14))
    ax0 = plt.subplot2grid((3, 2), (0, 0), projection='polar')
    ax1 = plt.subplot2grid((3, 2), (0, 1), projection='polar')
    file_ref = "Data/IceVolumeHudson_HE02.nc"
    time, vol = read_data(file_ref)
    surge_events = extract_surges(vol, thresh=-5e11)
    theta = 2 * np.pi / 1500.0 * time
    ax0.plot(theta[12:70], vol[12:70], color="black", lw=0.3)
    ax0.plot(theta[69:76], vol[69:76], color="black", lw=2)
    plot_do_cycle_polar(ax0)
    ax0.set_ylim(4e15, 5.9e15)
    ax0.set_theta_zero_location('N')
    ax0.set_theta_direction(-1)  # clockwise
    ax0.set_yticklabels([])
    ax0.yaxis.grid(False)
    ax0.set_xticks(xticks)
    ax0.tick_params(pad=5)
    ax0.set_xticklabels(labels, fontsize=13)
    ax0.xaxis.grid(True)
    create_subplot_labels(ax0, 0, 0)
    ax0.set_title("Single event cycle", fontsize=17, fontweight='bold')
    ax0.text(0.5, 0.80, 'warming phase', color='black',
             horizontalalignment='center', verticalalignment='center',
             transform=ax0.transAxes, fontsize=17)
    ax0.text(0.5, 0.20, 'cooling phase', color='black',
             horizontalalignment='center', verticalalignment='center',
             transform=ax0.transAxes, fontsize=17)
    ax0.text(-0.02, 0.20, 'Heinrich\nevent', color='blue',
             horizontalalignment='center', verticalalignment='center',
             transform=ax0.transAxes, fontsize=14)
    ax0.annotate('', xy=(0.13, 0.4), xycoords='axes fraction', xytext=(0, 0.25),
                 arrowprops=dict(arrowstyle='-|>', mutation_scale=25,
                                 color='blue', lw=2.5))
    ax0.text(1.12, 0.30, 'Ice stream\nbuild up', color='blue',
             horizontalalignment='center', verticalalignment='center',
             transform=ax0.transAxes, fontsize=14)
    ax0.annotate('', xy=(0.83, 0.3), xycoords='axes fraction',
                 xytext=(0.99, 0.3), arrowprops=dict(arrowstyle='-|>',
                                                     mutation_scale=25,
                                                     color='blue', lw=2.5))
    ax0.annotate('', xy=(0.68, 0.5), xycoords='axes fraction',
                 xytext=(1.05, 0.35), arrowprops=dict(arrowstyle='-|>',
                                                      mutation_scale=25,
                                                      color='blue', lw=2.5))
    ax1.plot(theta, vol, color="black", lw=0.3)
    for i in range(len(surge_events)):
        if i == 0:
            ax1.plot(theta[surge_events[i]], vol[surge_events[i]],
                     color='black', lw=2, label='Ice volume')
        else:
            ax1.plot(theta[surge_events[i]], vol[surge_events[i]],
                     color='black', lw=2)
    ax1.set_ylim(4e15, 5.9e15)
    ax1.set_theta_zero_location('N')
    ax1.set_theta_direction(-1)  # clockwise
    ax1.set_yticklabels([])
    ax1.yaxis.grid(False)
    ax1.set_xticks(xticks)
    ax1.tick_params(pad=5)
    ax1.set_xticklabels(labels, fontsize=13)
    ax1.legend(loc=1, bbox_to_anchor=(0.15, 0.1), fontsize=17)
    ax1.xaxis.grid(True)
    create_subplot_labels(ax1, 0, 1)
    ax1.set_title("Unforced reference simulation", fontsize=17, fontweight='bold')
    add_phase_locking_text(ax1, 1)
    count_row = 1
    count_col = 0

    for i_smb in range(len(runs[0])):
        for i_temp in range(len(runs)):
            if i_temp == 0 and i_smb == 0 or i_temp == 3 and i_smb == 0 or \
                    i_temp == 0 and i_smb == 5 or i_temp == 3 and i_smb == 5:

                fname = f'Data/IceVolume{args.region}_{runs[i_temp][i_smb]}.nc'
                print(fname)
                time, vol = read_data(fname)
                he_sample = get_he_events(runs[i_temp][i_smb], args.region, flux_thresh)
                sig_level = kuipers_test(he_sample, base_sample)
                plot_siglevel(sig_level, ax[count_row][count_col])
                surge_events = extract_surges(vol, thresh=grad_thresh)
                theta = 2 * np.pi / 1500.0 * time
                plot_he_lines(ax[count_row][count_col], theta, vol, surge_events)
                plot_do_cycle_polar(ax[count_row][count_col])
                set_axis_properties(args.region, ax[count_row][count_col],
                                    xticks, labels)
                create_subplot_labels(ax[count_row][count_col], count_row, count_col)
                add_phase_locking_text(ax[count_row][count_col], sig_level)
                if count_col == 0:
                    plot_panel_titles(ax[count_row][count_col],
                                      smb_pert[i_smb], 0)
                if count_row == 1:
                    plot_panel_titles(ax[count_row][count_col],
                                      t_pert[i_temp], 1)
                if count_col == 1:
                    count_col = 0
                    count_row += 1
                else:
                    count_col += 1
    plt.savefig(f'Fig02_Polar_{args.region}_Reduced.png', dpi=300,
                bbox_inches='tight')


if __name__ == '__main__':
    main()
