#!/usr/bin/env python
import sys
sys.path.append("../util")
from utility_fcts import *
import matplotlib.colors
from tueplots import figsizes


def main():
    # define some local vars for the plot
    label = ["a", "b", "c", "d"]
    snaps = [30, 63, 64, 84]
    inputfile = "Data/FlowLineData_Hudsoncoarse.nc"
    [thk, thk_grad, usurf, usurf_grad, bed, zb, temp_bottom, vel_base, tillw,
     bmelt, dbdt, time, taud, taub, distance] = read_flowline_vars(inputfile)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("",
                                                               ["lightgray",
                                                                "dimgray"])
    snowfall = ['medium', 'high', 'medium', 'small']

    plt.rcParams.update({"figure.dpi": 150})
    plt.rcParams.update(figsizes.icml2022_full(nrows=3, ncols=1))

    fig, ax = plt.subplots(nrows=4, ncols=1, sharex=True)

    for count, val in enumerate(snaps):
        zb_smooth = moving_average(zb[:, val])
        usurf_smooth = moving_average(usurf[:, val])
        bed_smooth = moving_average(bed[:, val])
        temp = temp_bottom[:, val]
        temp_layer = np.where((temp > -1e-02) * 1.0 == 1)
        start = temp_layer[0][0]
        end = temp_layer[0][-1]
        end = np.argwhere(np.isnan(temp))[0][0] - 3
        ocean_sindex = np.argmax(zb_smooth - bed_smooth > 2) - 1
        ax[count].fill_between(distance - 100, -500 * np.ones(len(distance)),
                               np.zeros(len(distance)), color='lightblue')
        polygon = ax[count].fill_between(distance - 100, zb_smooth, usurf_smooth,
                                         color='none')
        verts = np.vstack([p.vertices for p in polygon.get_paths()])
        gradient = ax[count].imshow(np.linspace(1, 0, 256).reshape(-1, 1),
                                    cmap=cmap, aspect='auto',
                                    extent=[verts[:, 0].min(),
                                            verts[:, 0].max(),
                                            verts[:, 1].min(),
                                            verts[:, 1].max()], zorder=1)
        gradient.set_clip_path(polygon.get_paths()[0], transform=ax[count].transData)
        ax[count].fill_between(distance[start - 1:end] - 100,
                               zb_smooth[start - 1:end] - 10,
                               zb_smooth[start - 1:end] + 100, color='red')
        ax[count].fill_between(distance[0:start] - 100, zb_smooth[0:start] - 10,
                               zb_smooth[0:start] + 100, color='blue')
        ax[count].fill_between(distance[ocean_sindex:] - 100,
                               zb_smooth[ocean_sindex:],
                               bed_smooth[ocean_sindex:], color='lightblue')
        ax[count].fill_between(distance - 100, bed_smooth,
                               -5000 * np.ones(len(distance)),
                               color='saddlebrown')
        ax[count].set_xlim([0, 2200])
        ax[count].set_ylim([-1000, 3500])
        ax[count].set_yticks([-1000, 0, 1000, 2000])
        ax[count].set_xticks([0, 1000, 2000])
        ax[count].set_xticklabels([])
        ax[count].set_ylabel("Elevation [m]")
        ax[count].text(0.015, 0.98, label[count], fontweight='bold',
                       transform=ax[count].transAxes, fontsize=12,
                       verticalalignment='top')
        if count > 0:
            ax[count].fill_between(distance - 100, zb_smooth_prev,
                                   usurf_smooth_prev, color='none',
                                   edgecolor='black', linestyle='dashed',
                                   linewidth=1.5)
        usurf_smooth_prev = usurf_smooth
        zb_smooth_prev = zb_smooth
        axins2 = inset_axes(ax[count], width="65%", height="52%",
                            bbox_to_anchor=(.41, .2, .55, .8),
                            bbox_transform=ax[count].transAxes, loc=1)
        create_do_forcing_plot(fig, axins2, count)
        draw_cloud(0.2, 0.95, ax[count])
        draw_cloud(0.4, 0.95, ax[count])
        draw_snow(0.2, 0.95, ax[count], snowfall[count])
        draw_snow(0.4, 0.95, ax[count], snowfall[count])
        if count == 3:
            ax[count].set_xticklabels([0, 1000, 2000])
            ax[count].set_xlabel("Distance [km]")
    plt.savefig("Fig01_Flowline.png", bbox_inches='tight', dpi=300)
    plt.show()
# sys.exit()


if __name__ == '__main__':
    main()
