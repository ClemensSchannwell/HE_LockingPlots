#! /usr/bin/env python

import argparse
import sys

sys.path.append("../util")
from utility_fcts import *


def create_subplot_labels(ax, row, col):
    label_mat = [
        ["(a)", "(b)", "(c)", "(d)"],
        ["(e)", "(f)", "(g)", "(h)"],
        ["(i)", "(j)", "(k)", "(l)"],
        ["(m)", "(n)", "(o)", "(p)"],
        ["(q)", "(r)", "(s)", "(t)"],
        ["(u)", "(v)", "(w)", "(x)"],
        ["(y)", "(z)", "(aa)", "(bb)"],
    ]
    ax.text(
        0.05,
        0.95,
        label_mat[row][col],
        color="black",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        fontsize=14,
    )

def plot_panel_titles_smb(ax, i_smb, i_temp, t_pert, smb_pert):
    if i_temp == 0:
        ax.set_title(
            "$\Delta$SMB = " + str(round(smb_pert[i_smb] / 1000, 3)) + " m/yr",
            fontsize=14,
            fontweight="bold",
        )
    if i_smb == 0:
        smb_str = (
            "$\Delta$Temp = " + str(t_pert[i_temp]) + "°C"
        )
        ax.text(
            -0.2,
            0.5,
            smb_str,
            color="black",
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            rotation=90,
            fontsize=14,
            fontweight="bold",
        )

def plot_panel_titles_temp(ax, i_smb, i_temp, t_pert, smb_pert):
    if i_smb == 0:
        ax.set_title(
            "$\Delta$Temp = " + str(t_pert[i_temp]) + "°C",
            fontsize=14,
            fontweight="bold",
        )
    if i_temp == 0:
        smb_str = (
            "$\Delta$SMB = " + str(round(smb_pert[i_smb] / 1000, 3)) + " m/yr"
        )
        ax.text(
            -0.2,
            0.5,
            smb_str,
            color="black",
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            rotation=90,
            fontsize=14,
            fontweight="bold",
        )


def set_axis_properties(region, ax, xticks, labels):
    if region == "Kenzie":
        ax.set_ylim(3.85e15, 4.95e15)
    else:
        ax.set_ylim(4e15, 5.9e15)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_yticklabels([])
    ax.yaxis.grid(False)
    ax.set_xticks(xticks)
    ax.set_xticklabels(labels)
    ax.xaxis.grid(True)
    ax.tick_params(pad=2)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-var",
        "--variable",
        choices=["temp", "smb"],
        help="Specify which plot you want to create",
        required=True,
    )
    args = parser.parse_args()
    clen = 1500
    grad_thresh = -3e11
    flux_thresh = 2.25e12
    labels = ["0", "300", "600", "900", "1200"]
    if args.variable == "temp":
        t_pert = [14, 10, 6]
        smb_pert = [100, 100, 100]
        runs = ["HE116", "HE117", "HE118"]
    else:
        t_pert = [0, 0, 0]
        smb_pert = [100, 75, 50]
        runs = ["HE119", "HE120", "HE121"]
    base_sample = get_he_events("HE02", "Hudson", clen=clen)
    ticks = np.array([0, 300, 600, 900, 1200])
    xticks = 2 * np.pi / clen * ticks
    m, n = 3, 1
    fig, ax = plt.subplots(
        n,
        m,
        subplot_kw={"projection": "polar"},
        gridspec_kw={
            "wspace": 0.43,
            "hspace": 0.175,
            "top": 0.99,
            "bottom": 0.01,
            "left": 0.03,
            "right": 0.97,
        },
        figsize=(12, 5),
    )
    for i_plot in range(len(runs)):
        fname = f"Data/IceVolumeHudson_{runs[i_plot]}.nc"
        time, vol = read_data(fname)
        he_sample = get_he_events(
            runs[i_plot], "Hudson", flux_thresh, clen=clen
        )
        sig_level = kuipers_test(he_sample, base_sample)
        plot_siglevel(sig_level, ax[i_plot])
        surge_events = extract_surges(vol, thresh=grad_thresh)
        theta = 2 * np.pi / clen * time
        plot_he_lines(ax[i_plot], theta, vol, surge_events)
        plot_do_cycle_polar(ax[i_plot], clen=clen)
        set_axis_properties(
            "Hudson", ax[i_plot], xticks, labels
        )
        if args.variable == "temp":
            plot_panel_titles_temp(ax[i_plot], 0, i_plot, t_pert, smb_pert)
        else:
            plot_panel_titles_smb(ax[i_plot], i_plot, 0, t_pert, smb_pert)
        create_subplot_labels(ax[i_plot], 0, i_plot)
    ax[0].legend(loc=1, bbox_to_anchor=(1.25, 1.1))
    plt.savefig(
        f"SIFig03_{args.variable}Only.png",
        dpi=300,
        bbox_inches="tight",
    )


if __name__ == "__main__":
    main()
