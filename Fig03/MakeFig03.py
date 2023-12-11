#! /usr/bin/env python
import numpy as np
import matplotlib.patches as mpatches
import sys
sys.path.append("../util")
from utility_fcts import *


def main():
    runs = [["HE68", "HE69", "HE70", "HE71", "HE72", "HE73"],
            ["HE46", "HE61", "HE55", "HE64", "HE49", "HE52"],
            ["HE47", "HE62", "HE56", "HE65", "HE50", "HE53"],
            ["HE48", "HE63", "HE57", "HE66", "HE51", "HE54"]]
    year_thresh_vec = [300, 600]
    sync_mat = compute_sync_mat(year_thresh_vec, runs)
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    axis_fs = 18
    tick_label_fs = 17
    colors = ["red", "black", "grey", "blue"]
    labels = ['300', '600']
    leg_labels = ["$\Delta$Temp = 14°C", "$\Delta$Temp = 10°C",
                 "$\Delta$Temp = 6°C","$\Delta$Temp = 2°C" ]
    xticks = []
    counter = 0
    for i in range(len(sync_mat[:, 0])):
        x = np.arange(counter, counter + 4)
        ax.bar(x, sync_mat[i, :], color=colors, width=1.0)
        xticks.append(counter + 1.5)
        counter = counter + 5.5
    # create patches for legend
    patches = [mpatches.Patch(color=colors[i], label=leg_labels[i]) for i in range(len(colors))]
    ax.set_xticks(xticks)
    ax.set_xticklabels(labels)
    # add legend to plot
    ax.legend(loc=9, bbox_to_anchor=(0.40, 0.99), handles=patches,
              prop={'size': axis_fs - 1}, labelspacing=1.5)
    ax.set_xlabel("Synchronicity window [years]", fontsize=axis_fs)
    ax.set_ylabel("Probability [-]", fontsize=axis_fs)
    ax.tick_params(axis='both', labelsize=tick_label_fs)
    ax.set_title("Synchronicity between Hudson and Mackenzie", fontsize=19,
                 fontweight='bold')
    plt.savefig("Fig03_Snychronicity.png", dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    main()
