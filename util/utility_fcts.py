#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot as plt, patches
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import sys


def create_do_forcing_plot(fig, ax, counter):
    inds = [30, 200, 400, 1100]
    t, ampl = create_forcing_fct(1500, 0.7)
    ax.scatter(t[inds[counter]], ampl[inds[counter]], s=80,
               facecolors='black', edgecolors='black', zorder=10)
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
    lc = LineCollection(segments, cmap='RdBu_r', norm=norm, zorder=1)
    # Set the values used for colormapping
    lc.set_array(dydx)
    lc.set_linewidth(5)
    plt.plot(t, np.zeros((len(t))), linestyle='dashed', color='gray')
    line = ax.add_collection(lc)
    ax.set_xlim(t.min(), t.max())
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.set_ylim(-1.2, 1.2)
    cax = inset_axes(ax, width="3%", height="100%", loc='lower left',
                     bbox_to_anchor=(1.01, 0., 1, 1), bbox_transform=ax.transAxes,
                     borderpad=0,)
    cbar = fig.colorbar(line, cax=cax)
    cbar.set_ticks([-0.0008, 0.0030])
    cbar.set_ticklabels(['cooling', 'warming'])
    cbar.ax.tick_params(size=0, rotation=90, labelsize=tick_label_size)


def create_forcing_fct(length, dist_factor):
    period = length / 2.0
    time = np.linspace(0, length, 1500)
    sine_wave = time * np.pi / period
    y = 1 / dist_factor * (np.arctan(dist_factor * np.sin(sine_wave) /
                                     (1 - dist_factor * np.cos(sine_wave))))
    y = 1 / y.max() * y
    return time, y


def draw_cloud(x, y, ax, color='darkgrey'):
    circle1 = patches.Circle((x, y), radius=0.05, color=color,
                             transform=ax.transAxes)
    circle2 = patches.Circle((x + 0.02, y - 0.05), radius=0.05, color=color,
                             transform=ax.transAxes)
    circle3 = patches.Circle((x, y - 0.07), radius=0.05, color=color,
                             transform=ax.transAxes)
    circle4 = patches.Circle((x - 0.02, y - 0.05), radius=0.05, color=color,
                             transform=ax.transAxes)
    ax.add_patch(circle1)
    ax.add_patch(circle2)
    ax.add_patch(circle3)
    ax.add_patch(circle4)


def draw_snow(x, y, ax, snowfall, color='silver'):
    if snowfall == 'medium':
        xpoints = [x + 0.04, x - 0.02, x + 0.02, x + 0.01]
        ypoints = [y - 0.14, y - 0.16, y - 0.185, y - 0.22]
    elif snowfall == 'high':
        xpoints = [x + 0.04, x - 0.02, x + 0.02, x + 0.01, x - 0.03]
        ypoints = [y - 0.14, y - 0.16, y - 0.185, y - 0.22, y - 0.21]
    else:
        xpoints = [x, x - 0.02]
        ypoints = [y - 0.19, y - 0.21]
    for i in range(len(xpoints)):
        ax.text(xpoints[i], ypoints[i], '*', color=color,
                horizontalalignment='center', verticalalignment='center',
                transform=ax.transAxes, fontsize=15, fontweight='bold')


def get_first_zeros_incols(var):
    indrow = (var == 0).argmax(axis=0)
    indcol = np.arange(0, len(var[0, :]), 1)
    return indrow, indcol


def mask_by_value(var2mask, value=0):
    var2mask[var2mask == value] = np.nan
    return var2mask


def mask_by_var(var2mask, maskvar, value=0):
    var2mask[maskvar == value] = np.nan
    return var2mask


def moving_average(vec, window_size=10):
    window = np.ones(window_size) / window_size
    return np.convolve(vec, window, mode='same')


def proc_flowline_vars(zs, zb, temp_bottom, thk_grad, usurf_grad, vel_base,
                       thk):
    indrow, indcol = get_first_zeros_incols(zs)
    zs = mask_by_value(zs)
    zs = set_ind_to_zero(zs, indrow, indcol)
    zb = set_ind_to_zero(zb, indrow, indcol)
    temp_bottom = mask_by_var(temp_bottom, thk)
    thk_grad = mask_by_var(thk_grad, thk)
    usurf_grad = mask_by_var(usurf_grad, thk)
    vel_base = mask_by_var(vel_base, -2e9)
    return zs, zb, temp_bottom, thk_grad, usurf_grad, vel_base


def read_flowline_vars(filepath):
    ncfile = Dataset(filepath)
    distance = ncfile.variables['y'][:]
    thk = ncfile.variables['thk'][:]
    thk_grad = np.array(np.gradient(thk, 10000, axis=0))
    usurf = ncfile.variables['usurf'][:]
    usurf_grad = np.array(np.gradient(usurf, 10000, axis=0))
    bed = ncfile.variables['topg'][:]
    zb = usurf - thk
    temp_bottom = ncfile.variables['temp_paBott'][:]
    vel_base = ncfile.variables['velbase_mag'][:]
    tillw = ncfile.variables['tillwat'][:]
    bmelt = ncfile.variables['bmelt'][:]
    dbdt = ncfile.variables['dbdt'][:]
    time = ncfile.variables['time'][:] / 86400 / 365
    taud = ncfile.variables['taud_mag'][:]
    taub = ncfile.variables['taub_mag'][:]
    [usurf, zb, temp_bottom, thk_grad, usurf_grad,vel_base] = proc_flowline_vars(usurf,zb, temp_bottom, thk_grad,
                                    usurf_grad, vel_base, thk)
    return (thk, thk_grad, usurf, usurf_grad, bed, zb, temp_bottom, vel_base,
            tillw, bmelt, dbdt, time, taud, taub, distance)


def set_ind_to_zero(var, indrow, indcol):
    var[indrow, indcol] = 0
    return var
