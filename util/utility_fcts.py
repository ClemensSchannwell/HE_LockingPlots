#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt, patches
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from astropy import stats as statsAs
from matplotlib.patches import Rectangle
import sys


def add_phase_locking_text(ax, sig_level):
    if sig_level >= 90:
        locking_str = "Phase locking"
        x = 0.85
        y = 1.00
    else:
        locking_str = "No phase locking"
        x = 0.9
        y = 1.0
    ax.text(
        x,
        y,
        locking_str,
        color="black",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        fontsize=13,
    )


def compute_flux(time, flux, threshold=2.25e12, clen=1500):
    clen_ind = int(clen / 100)
    cum_flux = np.zeros(clen_ind)
    cum_flux_max = np.zeros(clen_ind)
    he_timing = []
    counter = 0
    for i in range(0, len(time), clen_ind):
        if (i > clen_ind) and (i + clen_ind < len(time)) and counter == 0:
            flux_wind = flux[i : i + len(cum_flux)]
            time_wind = time[i : i + len(cum_flux)]
            # Con1=(FluxWind>Threshold)*1.0
            max_val = flux_wind.max()
            if max_val > threshold:
                counter = 1
                con2 = (flux_wind == max_val) * 1.0
                if (
                    np.nonzero(con2 == 1)[0][0] > 2
                    and np.nonzero(con2 == 1)[0][0] < 14
                ):
                    he_timing.append(time_wind[con2 == 1][0])
            else:
                con2 = np.zeros(clen_ind)
            cum_flux = cum_flux + (flux_wind > threshold) * 1.0
            cum_flux_max = cum_flux_max + con2
        else:
            counter = 0
    return cum_flux, cum_flux_max, he_timing


def compute_hist(time, flux, threshold=2.25e12, clen=1500):
    _, fmax, _ = compute_flux(time, flux, threshold=threshold, clen=clen)
    event = np.asarray(flux2hetiming(fmax))
    return event


def compute_sync_mat(thresh_vec, runs, clen=1500):
    counter = 0
    aggregate = np.array([])
    for thresh_count, thresh_val in enumerate(thresh_vec):
        sync_mat = np.zeros((len(runs[0][:]), len(runs[:])))
        for i_smb in range(len(runs[0])):
            for i_temp in range(len(runs)):
                fname_hud = f"Data/HudsonHov_{runs[i_temp][i_smb]}.nc"
                fname_ken = f"Data/KenzieHov_{runs[i_temp][i_smb]}.nc"
                thresh_hud = 2.25e12
                thresh_ken = 1.0e12
                time, flux = read_data(fname_hud, "flux_crossSection")
                _, _, he_timing_hud = compute_flux(
                    time, flux, threshold=thresh_hud, clen=clen
                )
                time, flux = read_data(fname_ken, "flux_crossSection")
                _, _, he_timing_ken = compute_flux(
                    time, flux, threshold=thresh_ken, clen=clen
                )
                event_counter = 0
                for he_count, he_val in enumerate(he_timing_hud):
                    diff = abs(he_val - he_timing_ken)
                    if np.sum(diff < thresh_val) > 0:
                        event_counter += 1

                sync_mat[i_smb, i_temp] = event_counter / len(he_timing_hud)
        if counter == 0:
            aggregate = np.mean(sync_mat, axis=0)
            counter = 1
        else:
            aggregate = np.vstack((aggregate, np.mean(sync_mat, axis=0)))
    return aggregate


def create_do_forcing_plot(fig, ax, counter):
    inds = [30, 200, 400, 1100]
    t, ampl = create_forcing_fct(1500, 0.7)
    ax.scatter(
        t[inds[counter]],
        ampl[inds[counter]],
        s=80,
        facecolors="black",
        edgecolors="black",
        zorder=10,
    )
    dydx = np.gradient(ampl, abs(t[0] - t[1]))
    # axis_fsize = 10
    tick_label_size = 9

    # Create a set of line segments so that we can color them individually
    # This creates the points as a N x 1 x 2 array so that we can stack points
    # together easily to get the segments. The segments array for line collection
    # needs to be (numlines) x (points per line) x 2 (for x and y)
    points = np.array([t, ampl]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    norm = plt.Normalize(-0.003, 0.003)
    lc = LineCollection(segments, cmap="RdBu_r", norm=norm, zorder=1)
    # Set the values used for colormapping
    lc.set_array(dydx)
    lc.set_linewidth(5)
    plt.plot(t, np.zeros((len(t))), linestyle="dashed", color="gray")
    line = ax.add_collection(lc)
    ax.set_xlim(t.min(), t.max())
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.set_ylim(-1.2, 1.2)
    cax = inset_axes(
        ax,
        width="3%",
        height="100%",
        loc="lower left",
        bbox_to_anchor=(1.01, 0.0, 1, 1),
        bbox_transform=ax.transAxes,
        borderpad=0,
    )
    cbar = fig.colorbar(line, cax=cax)
    cbar.set_ticks([-0.0008, 0.0030])
    cbar.set_ticklabels(["cooling", "warming"])
    cbar.ax.tick_params(size=0, rotation=90, labelsize=tick_label_size)


def create_forcing_fct(length, dist_factor):
    period = length / 2.0
    time = np.linspace(0, length, 1500)
    sine_wave = time * np.pi / period
    y = (
        1
        / dist_factor
        * (
            np.arctan(
                dist_factor
                * np.sin(sine_wave)
                / (1 - dist_factor * np.cos(sine_wave))
            )
        )
    )
    y = 1 / y.max() * y
    return time, y


def draw_cloud(x, y, ax, color="darkgrey"):
    circle1 = patches.Circle(
        (x, y), radius=0.05, color=color, transform=ax.transAxes
    )
    circle2 = patches.Circle(
        (x + 0.02, y - 0.05), radius=0.05, color=color, transform=ax.transAxes
    )
    circle3 = patches.Circle(
        (x, y - 0.07), radius=0.05, color=color, transform=ax.transAxes
    )
    circle4 = patches.Circle(
        (x - 0.02, y - 0.05), radius=0.05, color=color, transform=ax.transAxes
    )
    ax.add_patch(circle1)
    ax.add_patch(circle2)
    ax.add_patch(circle3)
    ax.add_patch(circle4)


def draw_snow(x, y, ax, snowfall, color="silver"):
    if snowfall == "medium":
        xpoints = [x + 0.04, x - 0.02, x + 0.02, x + 0.01]
        ypoints = [y - 0.14, y - 0.16, y - 0.185, y - 0.22]
    elif snowfall == "high":
        xpoints = [x + 0.04, x - 0.02, x + 0.02, x + 0.01, x - 0.03]
        ypoints = [y - 0.14, y - 0.16, y - 0.185, y - 0.22, y - 0.21]
    else:
        xpoints = [x, x - 0.02]
        ypoints = [y - 0.19, y - 0.21]
    for i in range(len(xpoints)):
        ax.text(
            xpoints[i],
            ypoints[i],
            "*",
            color=color,
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            fontsize=15,
            fontweight="bold",
        )


def extract_surges(ice_vol, time=100, thresh=-5e11):
    vol_grad = np.gradient(ice_vol, time)
    surge_inds = [i for i, n in enumerate(vol_grad) if n < thresh]
    surge_list = [[]]  # array of sub-arrays
    for i, num in enumerate(surge_inds):
        surge_list[-1].append(num)
        if (
            i != len(surge_inds) - 1
            and (surge_inds[i + 1] - surge_inds[i]) > 1
        ):
            surge_list.append([])
    filtered_surge_list = [x for x in surge_list if len(x) >= 4]
    return filtered_surge_list


def flux2hetiming(flux_vec, spacing=100):
    he_timing = []
    for i in range(len(flux_vec)):
        for j in range(int(flux_vec[i])):
            he_timing.append(i * spacing)
    return he_timing


def get_first_zeros_incols(var):
    indrow = (var == 0).argmax(axis=0)
    indcol = np.arange(0, len(var[0, :]), 1)
    return indrow, indcol


def get_he_events(run, region, threshold=2.25e12, clen=1500):
    fname = f"Data/{region}Hov_{run}.nc"
    time, flux = read_data(fname, "flux_crossSection")
    event = compute_hist(time, flux, clen=clen, threshold=threshold)
    return event


def kuipers_test(event, base_sample):
    k_test = abs(round(statsAs.kuiper_two(base_sample, event)[1], 3))
    return (1 - k_test) * 100


def mask_by_value(var2mask, value=0):
    var2mask[var2mask == value] = np.nan
    return var2mask


def mask_by_var(var2mask, maskvar, value=0):
    var2mask[maskvar == value] = np.nan
    return var2mask


def moving_average(vec, window_size=10):
    window = np.ones(window_size) / window_size
    return np.convolve(vec, window, mode="same")


def plot_do_cycle_polar(ax, clen=1500):
    n = 200
    if clen == 1000:
        theta_c1 = np.linspace(0, 2 * np.pi / clen * 100, n)
        theta_c2 = np.linspace(
            2 * np.pi / clen * 900, 2 * np.pi / clen * clen, n
        )
        theta_w1 = np.linspace(
            2 * np.pi / clen * 100, 2 * np.pi / clen * 500, n
        )
        theta_w2 = np.linspace(
            2 * np.pi / clen * 500, 2 * np.pi / clen * 900, n
        )
    elif clen == 2000:
        theta_c1 = np.linspace(0, 2 * np.pi / clen * 300, n)
        theta_c2 = np.linspace(
            2 * np.pi / clen * 1700, 2 * np.pi / clen * clen, n
        )
        theta_w1 = np.linspace(
            2 * np.pi / clen * 300, 2 * np.pi / clen * 1000, n
        )
        theta_w2 = np.linspace(
            2 * np.pi / clen * 1000, 2 * np.pi / clen * 1700, n
        )
    else:
        theta_c1 = np.linspace(0, 2 * np.pi / clen * 200, n)
        theta_c2 = np.linspace(
            2 * np.pi / clen * 1300, 2 * np.pi / 1500 * 1500, n
        )
        theta_w1 = np.linspace(
            2 * np.pi / clen * 200, 2 * np.pi / 1500 * 750, n
        )
        theta_w2 = np.linspace(
            2 * np.pi / clen * 750, 2 * np.pi / 1500 * 1300, n
        )
    t = np.linspace(0, 2 * np.pi, n)
    r = np.linspace(2.10e15, 5.9e15, 2)
    rg, tg = np.meshgrid(r, theta_c1)
    c = tg
    norm = mpl.colors.Normalize(0 - 1, theta_c1[-1])
    im = ax.pcolormesh(theta_c1, r, c.T, norm=norm, cmap="Reds_r")
    rg, tg = np.meshgrid(r, theta_w1)
    c = tg
    norm = mpl.colors.Normalize(theta_w1[0], theta_w1[-1] + 1)
    im = ax.pcolormesh(theta_w1, r, c.T, norm=norm, cmap="Blues")
    rg, tg = np.meshgrid(r, theta_w2)
    c = tg
    norm = mpl.colors.Normalize(theta_w2[0] - 1, theta_w2[-1])
    im = ax.pcolormesh(theta_w2, r, c.T, norm=norm, cmap="Blues_r")
    rg, tg = np.meshgrid(r, theta_c2)
    c = tg
    norm = mpl.colors.Normalize(theta_c2[0], theta_c2[-1] + 1)
    im = ax.pcolormesh(theta_c2, r, c.T, norm=norm, cmap="Reds")


def plot_he_lines(ax, theta, vol, surge_events):
    ax.plot(theta, vol, color="black", lw=0.3)
    for i in range(len(surge_events)):
        if i == 0:
            ax.plot(
                theta[surge_events[i]],
                vol[surge_events[i]],
                color="black",
                lw=2,
                label="Ice volume",
            )
        else:
            ax.plot(
                theta[surge_events[i]],
                vol[surge_events[i]],
                color="black",
                lw=2,
            )


def plot_panel_titles(ax, i_smb, i_temp):
    if i_temp == 0:
        smb_str = "$\Delta$SMB = " + str(round(i_smb / 1000, 3)) + " m/yr"
        ax.text(
            -0.2,
            0.5,
            smb_str,
            color="black",
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            rotation=90,
            fontsize=17,
            fontweight="bold",
        )
    elif i_temp == 1:
        ax.set_title(
            "$\Delta$Temp = " + str(i_smb) + "°C",
            fontsize=17,
            fontweight="bold",
        )


def plot_siglevel(sig_level, ax):
    alpha = 0.2
    rect = [0.1, 0.1, 0.8, 0.8]
    axx = ax.inset_axes(
        [0.0, 0.0, 1, 1],
    )
    autoAxis = axx.axis("off")
    if sig_level >= 99.0:
        rec = Rectangle(
            (-0.1, -0.05), 1.2, 1.1, fill=True, facecolor="lime", alpha=alpha
        )
        rec = axx.add_patch(rec)
        rec.set_clip_on(False)
    elif sig_level >= 90.0 and sig_level < 99.0:
        rec = Rectangle(
            (-0.1, -0.05), 1.2, 1.1, fill=True, facecolor="grey", alpha=alpha
        )
        rec = axx.add_patch(rec)
        rec.set_clip_on(False)
    axx.set_zorder(0)
    ax.set_zorder(axx.get_zorder() + 1)


def proc_flowline_vars(
    zs, zb, temp_bottom, thk_grad, usurf_grad, vel_base, thk
):
    indrow, indcol = get_first_zeros_incols(zs)
    zs = mask_by_value(zs)
    zs = set_ind_to_zero(zs, indrow, indcol)
    zb = set_ind_to_zero(zb, indrow, indcol)
    temp_bottom = mask_by_var(temp_bottom, thk)
    thk_grad = mask_by_var(thk_grad, thk)
    usurf_grad = mask_by_var(usurf_grad, thk)
    vel_base = mask_by_var(vel_base, -2e9)
    return zs, zb, temp_bottom, thk_grad, usurf_grad, vel_base


def read_data(fname, var="thk"):
    ncfile = Dataset(fname)
    time = ncfile.variables["time"][:] / (86400 * 365)
    time = time - time[0]
    if var == "thk":
        vol = ncfile.variables[var][:, 0, 0]
        return time, vol
    else:
        flux = ncfile.variables["flux_crossSection"][:]
        return time, flux


def read_flowline_vars(filepath):
    ncfile = Dataset(filepath)
    distance = ncfile.variables["y"][:]
    thk = ncfile.variables["thk"][:]
    thk_grad = np.array(np.gradient(thk, 10000, axis=0))
    usurf = ncfile.variables["usurf"][:]
    usurf_grad = np.array(np.gradient(usurf, 10000, axis=0))
    bed = ncfile.variables["topg"][:]
    zb = usurf - thk
    temp_bottom = ncfile.variables["temp_paBott"][:]
    vel_base = ncfile.variables["velbase_mag"][:]
    tillw = ncfile.variables["tillwat"][:]
    bmelt = ncfile.variables["bmelt"][:]
    dbdt = ncfile.variables["dbdt"][:]
    time = ncfile.variables["time"][:] / 86400 / 365
    taud = ncfile.variables["taud_mag"][:]
    taub = ncfile.variables["taub_mag"][:]
    [
        usurf,
        zb,
        temp_bottom,
        thk_grad,
        usurf_grad,
        vel_base,
    ] = proc_flowline_vars(
        usurf, zb, temp_bottom, thk_grad, usurf_grad, vel_base, thk
    )
    return (
        thk,
        thk_grad,
        usurf,
        usurf_grad,
        bed,
        zb,
        temp_bottom,
        vel_base,
        tillw,
        bmelt,
        dbdt,
        time,
        taud,
        taub,
        distance,
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
    ax.set_xticklabels(labels, fontsize=13)
    ax.xaxis.grid(True)
    ax.tick_params(pad=5)


def set_ind_to_zero(var, indrow, indcol):
    var[indrow, indcol] = 0
    return var
