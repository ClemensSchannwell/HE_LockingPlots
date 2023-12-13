#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as mpl
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
from matplotlib import ticker
import sys,os,glob
sys.path.append("/home/m/m300792/scripts")
import argparse
import cartopy.crs as ccrs
import cartopy
import cartopy.feature as cfeature
import colormaps as cmaps
import cmasher as cmr
# mpl.register_cmap(name='rignot',cmap=cmaps.rignot)
# mpl.register_cmap(name='rignot_r',cmap=cmaps.rignot_r)
# import colormaps as cmaps

def ReadNCVar(File,VarName):
    NCFile=Dataset(File)
    PlotVar = NCFile.variables[VarName][0,:,:]
    PlotVar[PlotVar==0]=np.nan
    Lat = NCFile.variables['lat'][:]
    Lon = NCFile.variables['lon'][:]
    Lon[Lon < 0] += 360
    X = NCFile.variables['x'][:]
    Y = NCFile.variables['y'][:]
    return PlotVar,Lat,Lon

File="/work/ba0989/m300792/HE_Runs/ExperimentsComposite/SMB_SensitivityExperiments/Template/pism_-049000/pism_-049000.nc"
FileDeltaT="/work/ba0989/m300792/HE_Runs/scripts/GenerateDOForcingT130/Delta_Strong_WeakHE.nc"
# FileDeltaT="DOState_WeakPISM.nc"
Var="thk"
Var1="ice_surface_temp"
Var2="climatic_mass_balance"
Thk, Lat, Lon = ReadNCVar(File,Var)
Temp, _, _ = ReadNCVar(FileDeltaT,Var1)
Temp = Temp.transpose()
SMB, _, _ = ReadNCVar(FileDeltaT,Var2)
SMB = SMB.transpose()

csLevs=[0,10,1000,2000,3000,4000,5000]
print(csLevs)
clevs=np.arange(0,4000,1000)
Tlevs=np.arange(-7,7,1)
cmapD = cmr.get_sub_cmap('seismic', 0.2, 0.8, N=9)
# cmap = plt.cm.get_cmap('seismic') # Get your favorite cmap
# new_cmap = sub_cmap(cmap, 0.0, 0.75)
# print(cmap)
# print(new_cmap)
# sys.exit()
# cmapD = plt.cm.get_cmap(new_cmap,9)
# print(clevs)

# fig,_ = plt.subplots(figsize=(8,4))
fig, ax = plt.subplots(1,2,figsize=(10,4),constrained_layout=True, subplot_kw={'projection': ccrs.NorthPolarStereo(central_longitude=-45)})
# ax[0] = plt.subplot(projection=ccrs.NorthPolarStereo(central_longitude=-45))
ax[0].coastlines(color="grey",zorder=1,resolution='110m')
ax[0].set_extent([-120, 30, 50, 90], crs=ccrs.PlateCarree())
ax[0].gridlines(draw_labels=True, linewidth=0.5, linestyle='--', color='black')
pc=ax[0].pcolormesh(Lon, Lat, SMB*3.155e4,transform=ccrs.PlateCarree(), cmap=cmapD,
            vmin=-0.1125,vmax=0.1125)
ax[0].contour(Lon,Lat,Thk,clevs,transform=ccrs.PlateCarree(),colors=['black'])
cbar = fig.colorbar(pc,ax=ax[0],fraction=0.036,pad=0.02)
cbar.set_label('$\Delta$SMB [m/yr]')
# plt.show()
# sys.exit()
# plt.savefig("DeltaSMBT130HE.png",bbox_inches='tight',dpi=300)

cmapD = cmr.get_sub_cmap('seismic', 0.2, 0.8, N=7)
# ax = plt.subplot(projection=ccrs.NorthPolarStereo(central_longitude=-45))
ax[1].coastlines(color="grey",zorder=1,resolution='110m')
ax[1].set_extent([-120, 30, 50, 90], crs=ccrs.PlateCarree())
ax[1].gridlines(draw_labels=True, linewidth=0.5, linestyle='--', color='black')
pc=ax[1].pcolormesh(Lon, Lat, Temp,transform=ccrs.PlateCarree(), cmap=cmapD,
            vmin=-7,vmax=7)
ax[1].contour(Lon,Lat,Thk,clevs,transform=ccrs.PlateCarree(),colors=['black'])
cbar = fig.colorbar(pc,ax=ax[1],fraction=0.036,pad=0.02)
cbar.set_label('$\Delta$T [K]')
plt.savefig("DeltaTempT130HE.png",bbox_inches='tight',dpi=300)
plt.show()
sys.exit()
# Temp = PlotMask(MaskFiles[i],Colors[i],Lab[i])
    # print(Temp.collections[0])
    # plt.show()
    # sys.exit()
    # CS[i] = Temp.collections[0]
# plt.legend(CS,Lab)
# plt.show()
# print(CS)
# sys.exit()
# # ax.contour(Lon,Lat,Mask,colors=['red'],transform=ccrs.PlateCarree())
# cbar = fig.colorbar(pc)
# # plt.savefig("CartopyMap.png",bbox_inches='tight',dpi=300)
# plt.show()
# sys.exit()
# print(X)
# ax.contourf(Lon,Lat,PlotVar,clevs,cmap='Blues',transform=ccrs.PlateCarree())
