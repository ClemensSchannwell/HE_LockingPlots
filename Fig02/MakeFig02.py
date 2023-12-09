#! /usr/bin/env python

import argparse
import sys
sys.path.append("../util")
from utility_fcts import *


def create_subplot_labels(ax, row, col):
    label_mat = [['(a)', '(b)'], ['(c)', '(d)'], ['(e)', '(f)']]
    ax.text(0.05, 0.95, label_mat[row][col], color='black',
            horizontalalignment='center', verticalalignment='center',
            transform=ax.transAxes, fontsize=17)

# import numpy as np
# import numpy.random
# import copy
# import matplotlib.pyplot as plt
# import scipy.stats as st
# from scipy.stats import bootstrap
# from scipy import integrate
# from scipy import stats
# from astropy import stats as statsAs
# from netCDF4 import Dataset
# import sys
# import matplotlib as mpl
# from matplotlib.patches import Rectangle


# def add_phase_locking_text(ax, sig_level):
    # print(sig_level)
    # if sig_level >= 90:
        # locking_str = "Phase locking"
        # x = 0.85
        # y = 1.00
    # else:
        # locking_str = "No phase locking"
        # x = 0.9
        # y = 1.0
    # ax.text(x,y, locking_str ,color='black',
            # horizontalalignment='center', verticalalignment='center',
            # transform=ax.transAxes,fontsize=13)


# def compute_flux(Time,Flux,Threshold=2.25e12,Clen=1500):
    # CLenInd=int(Clen/100)
    # # print(Threshold)
    # CumFlux=np.zeros(CLenInd)
    # CumFluxMax=np.zeros(CLenInd)
    # counter=0
    # for i in range(0,len(Time),CLenInd):
        # # print('OuterLoop:',Time[i],Time[i+len(CumFlux)])
        # if (i > CLenInd) and (i+CLenInd < len(Time)) and counter == 0:
            # FluxWind=Flux[i:i+len(CumFlux)]
            # Con1=(FluxWind>Threshold)*1.0
            # MaxVal=FluxWind.max()
            # # print('Inner loop:',Time[i],Time[i+len(CumFlux)])
            # # sys.exit()
            # if MaxVal > Threshold:
                # # print(i,i+len(CumFlux))
                # # print(Time[i]-65000,Time[i+len(CumFlux)]-65000)
                # # print(FluxWind)
                # # sys.exit()
                # counter=1
                # Con2=(FluxWind==MaxVal)*1.0
            # else:
                # Con2=np.zeros(CLenInd)
            # CumFlux=CumFlux+ (FluxWind>Threshold)*1.0
            # CumFluxMax=CumFluxMax + Con2
        # else:
            # counter = 0
    # return CumFlux,CumFluxMax

# def compute_hist(Time,Flux,Threshold=2.25e12,CLen=1500):
    # _,FMax = compute_flux(Time,Flux,Threshold=Threshold,Clen=CLen)
    # # print(FMax)
    # Event = np.asarray(flux2hetiming(FMax))
    # return Event

# def get_he_events(Run,Threshold=2.25e12,CLen=1500):
    # File=BasePath + Run + "/Postprocessing/" + Region + "Hov.nc"
    # T, F = ReadData(File,'flux_crossSection')
    # Event = compute_hist(T,F,CLen=CLen, Threshold=Threshold)
    # return Event

# def extract_surges(icevol,time=100,thresh=-5e11):
    # vol_grad = np.gradient(icevol,time)
    # # fig,ax1 = plt.subplots(1,1)
    # # print(vol_grad)
    # # ax1.plot(T,icevol)
    # surge_inds = [i for i, n in enumerate(vol_grad) if n < thresh]
    # surge_list = [[]] # array of sub-arrays
    # for i, num in enumerate(surge_inds):          # go through each element after the first
        # surge_list[-1].append(num)                  # Add it to the last sub-array
        # if i != len(surge_inds)-1 and (surge_inds[i+1] - surge_inds[i]) > 1 :   # if 5 encountered and not last element
            # surge_list.append([])
    # filtered_surge_list=[x for x in surge_list if len(x)>=4]
# # print(filtered_surge_list)
# # for j in range(len(filtered_surge_list)):
    # # ax1.plot(T[filtered_surge_list[j]],icevol[filtered_surge_list[j]], color='green', lw=3)
# # plt.show()
# # sys.exit()
# return filtered_surge_list

# def flux2hetiming(FluxVec,spacing=100):
# HETiming = []
# # print(FluxVec)
# for i in range(len(FluxVec)):
    # for j in range(int(FluxVec[i])):
        # HETiming.append(i*spacing)
# return HETiming

# def kuipers_test(Event,BaseSample,FS=8):
# KTest = abs(round(statsAs.kuiper_two(BaseSample,Event)[1],3))
# KTestStr = "p-value: " + str(KTest)
# return (1-KTest)*100

# def plot_do_cycle_polar(ax):
# n = 200
# thetaC1 = np.linspace(0, 2*np.pi / 1500*200, n)
# thetaC2 = np.linspace(2*np.pi / 1500*1300,2*np.pi/1500*1500, n)
# thetaW1 = np.linspace(2*np.pi / 1500*200,2*np.pi/1500*750, n)
# thetaW2 = np.linspace(2*np.pi / 1500*750,2*np.pi/1500*1300, n)
# t = np.linspace(0,2*np.pi,n)   #theta values
# r = np.linspace(2.10e15,5.9e15,2)        #radius values full circle
# rg, tg = np.meshgrid(r,thetaC1)      #create a r,theta meshgrid
# c = tg                         #define color values as theta value
# norm = mpl.colors.Normalize(0-1, thetaC1[-1])
# im = ax.pcolormesh(thetaC1, r, c.T, norm=norm, cmap='Reds_r')
# rg, tg = np.meshgrid(r,thetaW1)
# c = tg
# norm = mpl.colors.Normalize(thetaW1[0], thetaW1[-1]+1)
# im = ax.pcolormesh(thetaW1, r, c.T,norm=norm, cmap='Blues')
# rg, tg = np.meshgrid(r,thetaW2)      #create a r,theta meshgrid
# c = tg                         #define color values as theta value
# norm = mpl.colors.Normalize(thetaW2[0]-1, thetaW2[-1])
# im = ax.pcolormesh(thetaW2, r, c.T,norm=norm, cmap='Blues_r')
# rg, tg = np.meshgrid(r,thetaC2)      #create a r,theta meshgrid
# c = tg                         #define color values as theta value
# norm = mpl.colors.Normalize(thetaC2[0], thetaC2[-1]+1)
# im = ax.pcolormesh(thetaC2, r, c.T,norm=norm, cmap='Reds')

# def plot_he_lines(ax, Theta, Vol, surge_events):
# ax.plot(Theta,Vol,color="black",lw=0.3)
# for i in range(len(surge_events)):
    # if i == 0:
        # ax.plot(Theta[surge_events[i]],Vol[surge_events[i]], color='black', lw=2,
            # label='Ice volume')
    # else:
        # ax.plot(Theta[surge_events[i]],Vol[surge_events[i]], color='black', lw=2)

# def plot_panel_titles(ax,iSMB,iTemp):
    # # ax[iSMB][iTemp].set_title("$\Delta$Temp = " + str(TPert[iTemp])+"°C", \
        # # fontsize=14,fontweight='bold')
# if iTemp == 0:
    # SMBStr = "$\Delta$SMB = " + str(round(iSMB/1000,3))+" m/yr"
    # ax.text(-0.2, 0.5, SMBStr,color="black",
    # horizontalalignment='center', verticalalignment='center',
    # transform=ax.transAxes, rotation=90,
    # fontsize=17,fontweight='bold')
# elif iTemp == 1:
    # ax.set_title("$\Delta$Temp = " + str(iSMB)+"°C", \
        # fontsize=17,fontweight='bold')

# def plot_siglevel(SigLevel,ax):
# Alpha=0.2
# rect  = [0.1,0.1,0.8,0.8]
# axx = ax.inset_axes([.0, .0, 1, 1],)
# autoAxis = axx.axis('off')
# if SigLevel >= 99.0:
    # rec = Rectangle((-0.1,-0.05),1.2,1.1,fill=True,
            # facecolor='lime',alpha=Alpha)
    # rec = axx.add_patch(rec)
    # rec.set_clip_on(False)
# elif SigLevel >= 90.0 and SigLevel < 99.0:
    # rec = Rectangle((-0.1,-0.05),1.2,1.1,fill=True,
            # facecolor='grey',alpha=Alpha)
    # rec = axx.add_patch(rec)
    # rec.set_clip_on(False)
# axx.set_zorder(0)
# ax.set_zorder(axx.get_zorder()+1)

# def ReadData(FilePath,var='thk'):
# NCFile=Dataset(FilePath)
# Time=NCFile.variables['time'][:]/(86400*365)
# Time = Time - Time[0]
# if var == 'thk':
    # Vol=NCFile.variables[var][:,0,0]
    # return Time,Vol
# else:
    # Flux=NCFile.variables['flux_crossSection'][:]
    # return Time,Flux
# # [print(FilePath,Time[x]) for x in range(len(Time))]

# def set_axis_properties(ax,xticks,labels):
# if Region == "Kenzie":
    # ax.set_ylim(3.85e15,4.95e15)
# else:
    # ax.set_ylim(4e15,5.9e15)
# ax.set_theta_zero_location('N')
# ax.set_theta_direction(-1) # clockwise
# ax.set_yticklabels([])
# ax.yaxis.grid(False)
# ax.set_xticks(xticks)
# ax.set_xticklabels(labels,fontsize=13)
# ax.xaxis.grid(True)
# ax.tick_params(pad=5)

def main():
    # global base_path, region, t_pert, smb_pert
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
    t_pert = [14, 10, 6, 2]
    smb_pert = [150, 100, 75, 50, 37.5, 9.375]
    runs = [["HE68", "HE69", "HE70", "HE71", "HE72", "HE73"],
            ["HE46", "HE61", "HE55", "HE64", "HE49", "HE52"],
            ["HE47", "HE62", "HE56", "HE65", "HE50", "HE53"],
            ["HE48", "HE63", "HE57", "HE66", "HE51", "HE54"]]
    ticks = np.array([0, 300, 600, 900, 1200])
    xticks = 2 * np.pi / 1500.0 * ticks
    labels = ['0', '300', '600', '900', '1200']
    m, n = 3, 2
    fig, ax = plt.subplots(m, n, subplot_kw={'projection': 'polar'},
                           gridspec_kw={'wspace': 0.15, 'hspace': 0.175,
                                        'top': 0.99, 'bottom': 0.01,
                                        'left': 0.03, 'right': 0.97},
                           figsize=(10.5, 14))
    ax0 = plt.subplot2grid((3, 2), (0, 0), projection='polar')
    ax1 = plt.subplot2grid((3, 2), (0, 1), projection='polar')
    file_ref = "Data/IceVolumeHudson_HE02.nc"
    time, vol = read_data(file_ref)
    surge_events = extract_surges(vol, thresh=-5e11)
    theta = 2 * np.pi / 1500.0 * time
    ax0.plot(theta[12:70], vol[12:70], color="black", lw=0.3)
    ax0.plot(theta[69:76], vol[69:76], color="black", lw=2)
    plot_do_cycle_polar(ax0)
    ax0.set_ylim(4e15, 5.9e15)
    ax0.set_theta_zero_location('N')
    ax0.set_theta_direction(-1)  # clockwise
    ax0.set_yticklabels([])
    ax0.yaxis.grid(False)
    ax0.set_xticks(xticks)
    ax0.tick_params(pad=5)
    ax0.set_xticklabels(labels, fontsize=13)
    ax0.xaxis.grid(True)
    create_subplot_labels(ax0, 0, 0)
    ax0.set_title("Single event cycle", fontsize=17, fontweight='bold')
    ax0.text(0.5, 0.80, 'interstadial phase', color='black',
             horizontalalignment='center', verticalalignment='center',
             transform=ax0.transAxes, fontsize=17)
    ax0.text(0.5, 0.20, 'stadial phase', color='black',
             horizontalalignment='center', verticalalignment='center',
             transform=ax0.transAxes, fontsize=17)
    ax0.text(-0.02, 0.20, 'Heinrich\nevent', color='blue',
             horizontalalignment='center', verticalalignment='center',
             transform=ax0.transAxes, fontsize=14)
    ax0.annotate('', xy=(0.13, 0.4), xycoords='axes fraction', xytext=(0, 0.25),
                 arrowprops=dict(arrowstyle='-|>', mutation_scale=25,
                                 color='blue', lw=2.5))
    ax0.text(1.12, 0.30, 'Ice stream\nbuild up', color='blue',
             horizontalalignment='center', verticalalignment='center',
             transform=ax0.transAxes, fontsize=14)
    ax0.annotate('', xy=(0.83, 0.3), xycoords='axes fraction',
                 xytext=(0.99, 0.3), arrowprops=dict(arrowstyle='-|>',
                                                     mutation_scale=25,
                                                     color='blue', lw=2.5))
    ax0.annotate('', xy=(0.68, 0.5), xycoords='axes fraction',
                 xytext=(1.05, 0.35), arrowprops=dict(arrowstyle='-|>',
                                                      mutation_scale=25,
                                                      color='blue', lw=2.5))
    ax1.plot(theta, vol, color="black", lw=0.3)
    for i in range(len(surge_events)):
        if i == 0:
            ax1.plot(theta[surge_events[i]], vol[surge_events[i]],
                     color='black', lw=2, label='Ice volume')
        else:
            ax1.plot(theta[surge_events[i]], vol[surge_events[i]],
                     color='black', lw=2)
    ax1.set_ylim(4e15, 5.9e15)
    ax1.set_theta_zero_location('N')
    ax1.set_theta_direction(-1)  # clockwise
    ax1.set_yticklabels([])
    ax1.yaxis.grid(False)
    ax1.set_xticks(xticks)
    ax1.tick_params(pad=5)
    ax1.set_xticklabels(labels, fontsize=13)
    ax1.legend(loc=1, bbox_to_anchor=(0.15, 0.1), fontsize=17)
    ax1.xaxis.grid(True)
    create_subplot_labels(ax1, 0, 1)
    ax1.set_title("Unforced reference simulation", fontsize=17, fontweight='bold')
    add_phase_locking_text(ax1, 1)
    count_row = 1
    count_col = 0

    for i_smb in range(len(runs[0])):
        for i_temp in range(len(runs)):
            if i_temp == 0 and i_smb == 0 or i_temp == 3 and i_smb == 0 or \
                    i_temp == 0 and i_smb == 5 or i_temp == 3 and i_smb == 5:

                fname = f'Data/IceVolume{args.region}_{runs[i_temp][i_smb]}.nc'
                print(fname)
                time, vol = read_data(fname)
                he_sample = get_he_events(runs[i_temp][i_smb], args.region, flux_thresh)
                sig_level = kuipers_test(he_sample, base_sample)
                plot_siglevel(sig_level, ax[count_row][count_col])
                surge_events = extract_surges(vol, thresh=grad_thresh)
                theta = 2 * np.pi / 1500.0 * time
                plot_he_lines(ax[count_row][count_col], theta, vol, surge_events)
                plot_do_cycle_polar(ax[count_row][count_col])
                set_axis_properties(args.region, ax[count_row][count_col],
                                    xticks, labels)
                create_subplot_labels(ax[count_row][count_col], count_row, count_col)
                add_phase_locking_text(ax[count_row][count_col], sig_level)
                if count_col == 0:
                    plot_panel_titles(ax[count_row][count_col],
                                      smb_pert[i_smb], 0)
                if count_row == 1:
                    plot_panel_titles(ax[count_row][count_col],
                                      t_pert[i_temp], 1)
                if count_col == 1:
                    count_col = 0
                    count_row += 1
                else:
                    count_col += 1
    plt.savefig(f'Fig02_Polar_{args.region}_Reduced.png', dpi=300,
                bbox_inches='tight')


if __name__ == '__main__':
    main()
