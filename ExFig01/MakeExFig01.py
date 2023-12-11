#! /usr/bin/env python

import numpy as np
import sys
sys.path.append("../util")
from utility_fcts import *


def create_subplot_labels(ax, row):
    label_vec = ['(a)', '(b)']
    ax.text(-0.05, 1.05, label_vec[row], color='black',
            horizontalalignment='center', verticalalignment='center',
            transform=ax.transAxes, fontsize=17)


def main():
    ticks = np.array([0, 300, 600, 900, 1200])
    xticks = 2 * np.pi / 1500.0 * ticks
    labels = ['0', '300', '600', '900', '1200']
    fig, ax = plt.subplots(nrows=1, ncols=2, subplot_kw={'projection':
                                                         'polar'},
                           figsize=(11, 3))
    regions = ["Hudson", "Kenzie"]
    grad_thresh = [-3e11, -2e11]
    for count, ireg in enumerate(regions):
        fname = f"Data/IceVolume{ireg}_HE02.nc"
        time, vol = read_data(fname)
        surge_events = extract_surges(vol, thresh=grad_thresh[count])
        theta = 2 * np.pi / 1500.0 * time
        plot_he_lines(ax[count], theta, vol, surge_events)
        plot_do_cycle_polar(ax[count])
        set_axis_properties(ireg, ax[count], xticks, labels)
        create_subplot_labels(ax[count], count)
    ax[0].legend(loc=1, bbox_to_anchor=(2.02, 0.5), fontsize=17)
    plt.savefig('ExFig01_base_runs.png', dpi=300, bbox_inches='tight')
    plt.show()
    sys.exit()


if __name__ == '__main__':
    main()
