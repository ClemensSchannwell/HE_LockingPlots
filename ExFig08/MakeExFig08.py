#! /usr/bin/env python

import argparse
import sys
import numpy as np
sys.path.append("../util")
from utility_fcts import *


def create_subplot_labels(ax, row, col):
    label_mat = [['(a)', '(b)', '(c)', '(d)'], ['(e)', '(f)', '(g)', '(h)'],
                 ['(i)', '(j)', '(k)', '(l)'], ['(m)', '(n)', '(o)', '(p)'],
                 ['(q)', '(r)', '(s)', '(t)'], ['(u)', '(v)', '(w)', '(x)'],
                 ['(y)', '(z)', '(aa)', '(bb)']]
    ax.text(0.00, 1.00, label_mat[row][col], color='black',
            horizontalalignment='center', verticalalignment='center',
            transform=ax.transAxes, fontsize=14)


def get_he_timing(region, run, phase, surge_events):
    fname = f'Data/{region}Hov_{run}.nc'
    time, flux = read_data(fname, "flux_crossSection")
    max_ind = [surge_events[isurge][np.argmax(flux[surge_events[isurge]])]
               for isurge in range(len(surge_events))]
    return phase[max_ind]


def plot_panel_titles(ax, i_smb, i_temp, t_pert, smb_pert):
    if i_smb == 0:
        ax.set_title("$\Delta$Temp = " + str(t_pert[i_temp]) + "°C",
                     fontsize=14, fontweight='bold')
    if i_temp == 0:
        smb_str = "$\Delta$SMB = " + str(round(smb_pert / 1000, 3)) + " m/yr"
        ax.text(-0.2, 0.5, smb_str, color="black", horizontalalignment='center',
                verticalalignment='center', transform=ax.transAxes, rotation=90,
                fontsize=14, fontweight='bold')


def read_random_forcing():
    nc_file = Dataset("Data/random_forcing.nc")
    phase = nc_file.variables['phase'][0:887]
    theta = 2 * np.pi / 360 * phase
    return phase, theta


def set_axis_properties(region, ax):
    if region == "Kenzie":
        ax.set_ylim(3.85e15, 4.95e15)
    else:
        ax.set_ylim(4e15, 5.9e15)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_yticklabels([])
    ax.yaxis.grid(False)
    ax.xaxis.grid(True)
    ax.tick_params(pad=2)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-reg", "--region", choices=["Hudson", "Kenzie"],
                        help="Specify which plot you want to create",
                        required=True)

    t_pert = [14, 10, 6, 2]
    smb_pert = 100
    args = parser.parse_args()
    if args.region == "Kenzie":
        base_sample = get_he_events("HE02", args.region, threshold=1.0e12)
        grad_thresh = -2e11
        flux_thresh = 1.5e12
    else:
        base_sample = get_he_events("HE02", args.region)
        grad_thresh = -3e11
        flux_thresh = 2.25e12
    runs = ["HE114", "HE113", "HE111", "HE112"]
    xticks = np.array([0, 90, 135, 180, 225, 270, 315])
    labels = ['0', '90', '135', '180', '225', '270', '315']
    m, n = 4, 1
    fig, ax = plt.subplots(n, m, subplot_kw={'projection': 'polar'},
                           gridspec_kw={'wspace': 0.43, 'hspace': 0.175,
                                        'top': 0.99, 'bottom': 0.01,
                                        'left': 0.03, 'right': 0.97},
                           figsize=(12, 5))
    for i_mult in range(len(runs)):
        fname = f'Data/IceVolume{args.region}_{runs[i_mult]}.nc'
        time, vol = read_data(fname)
        surge_events = extract_surges(vol, thresh=grad_thresh)
        phase, theta = read_random_forcing()
        base_sample = 360 * base_sample / 1500
        plot_he_lines(ax[i_mult], theta, vol, surge_events)
        he_sample = get_he_timing(args.region, runs[i_mult], phase, surge_events)
        sig_level = kuipers_test(he_sample, base_sample)
        plot_siglevel(sig_level, ax[i_mult])
        plot_do_cycle_polar(ax[i_mult])
        set_axis_properties(args.region, ax[i_mult])
        plot_panel_titles(ax[i_mult], 0,  i_mult, t_pert, smb_pert)
        create_subplot_labels(ax[i_mult], 0, i_mult)
    ax[0].legend(loc=1, bbox_to_anchor=(1.37, 1.15))
    plt.savefig(f'ExFig08_Polar_{args.region}.png', dpi=300,
                bbox_inches='tight')


if __name__ == '__main__':
    main()
