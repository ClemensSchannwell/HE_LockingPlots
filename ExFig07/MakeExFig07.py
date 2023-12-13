#! /usr/bin/env python

import argparse
import sys
from tueplots import figsizes
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


def plot_panel_titles(ax, i_mult):
    forc_str = f"x{i_mult + 1}"
    ax.text(1.05, 0.90, forc_str, color="black", horizontalalignment='center',
            verticalalignment='center', transform=ax.transAxes, fontsize=14)


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
    else:
        base_sample = get_he_events("HE02", args.region)
        grad_thresh = -3e11
        flux_thresh = 2.25e12
        GradThresh = -2e11
    runs = ["HE97", "HE90", "HE91", "HE92"]
    ticks = np.array([0, 300, 600, 900, 1200])
    xticks = 2 * np.pi / 1500.0 * ticks
    labels = ['0', '300', '600', '900', '1200']
    plt.rcParams.update({"figure.dpi": 150})
    plt.rcParams.update(figsizes.icml2022_half(nrows=len(runs), ncols=1))
    fig, ax = plt.subplots(len(runs), 1, subplot_kw={'projection': 'polar'})

    for i_mult in range(len(runs)):
        fname = f'Data/IceVolume{args.region}_{runs[i_mult]}.nc'
        time, vol = read_data(fname)
        he_sample = get_he_events(runs[i_mult], args.region, flux_thresh)
        sig_level = kuipers_test(he_sample, base_sample)
        plot_siglevel(sig_level, ax[i_mult])
        surge_events = extract_surges(vol, thresh=grad_thresh)
        theta = 2 * np.pi / 1500.0 * time
        plot_he_lines(ax[i_mult], theta, vol, surge_events)
        plot_do_cycle_polar(ax[i_mult])
        set_axis_properties(args.region, ax[i_mult], xticks, labels)
        plot_panel_titles(ax[i_mult], i_mult)
        create_subplot_labels(ax[i_mult], 0, i_mult)
    ax[0].legend(loc=1, bbox_to_anchor=(1.4, 1.2))
    plt.savefig(f'ExFig07_Polar_{args.region}.png', dpi=300,
                bbox_inches='tight')


if __name__ == '__main__':
    main()
