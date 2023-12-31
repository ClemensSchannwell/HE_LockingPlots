#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cmasher as cmr

# import matplotlib.pylab as mpl
# import matplotlib.colors as colors
# from matplotlib.colors import LogNorm
# from matplotlib import ticker
# import sys,os,glob
# sys.path.append("/home/m/m300792/scripts")
# import argparse
# import cartopy.crs as ccrs
# import cartopy
# import cartopy.feature as cfeature
# import colormaps as cmaps
# import cmasher as cmr


def read_var(fname, var_name):
    nc_file = Dataset(fname)
    pl_var = nc_file.variables[var_name][0, :, :]
    pl_var[pl_var == 0] = np.nan
    lat = nc_file.variables["lat"][:]
    lon = nc_file.variables["lon"][:]
    lon[lon < 0] += 360
    # X = NCFile.variables['x'][:]
    # Y = NCFile.variables['y'][:]
    return pl_var, lat, lon


def main():
    fname0 = "Data/pism_thk.nc"
    thk, lat, lon = read_var(fname0, "thk")
    fname1 = "Data/Delta_Strong_WeakHE.nc"
    (
        temp,
        _,
        _,
    ) = read_var(fname1, "ice_surface_temp")
    temp = temp.transpose()
    (
        smb,
        _,
        _,
    ) = read_var(fname1, "climatic_mass_balance")
    smb = smb.transpose() * 3.1554e4

    fig, ax = plt.subplots(
        1,
        2,
        figsize=(10, 4),
        constrained_layout=True,
        subplot_kw={
            "projection": ccrs.NorthPolarStereo(central_longitude=-45)
        },
    )
    clevs = np.arange(0, 4000, 1000)
    pl_vars = [smb, temp]
    var_max = [0.1125, 7]
    cm_classes = [9, 7]
    labels = ["$\Delta$SMB [m/yr]", "$\Delta$T [K]"]

    for i_count, i_ax in enumerate(ax):
        cmapD = cmr.get_sub_cmap("seismic", 0.2, 0.8, N=cm_classes[i_count])
        i_ax.coastlines(color="grey", zorder=1, resolution="110m")
        i_ax.set_extent([-120, 30, 50, 90], crs=ccrs.PlateCarree())
        i_ax.gridlines(
            draw_labels=True, linewidth=0.5, linestyle="--", color="black"
        )
        pc = i_ax.pcolormesh(
            lon,
            lat,
            pl_vars[i_count],
            transform=ccrs.PlateCarree(),
            cmap=cmapD,
            vmin=-var_max[i_count],
            vmax=var_max[i_count],
        )
        i_ax.contour(
            lon,
            lat,
            thk,
            clevs,
            transform=ccrs.PlateCarree(),
            colors=["black"],
        )
        cbar = fig.colorbar(pc, ax=i_ax, fraction=0.036, pad=0.02)
        cbar.set_label(labels[i_count])
    plt.savefig(f"ExFig06_DO2D.png", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
