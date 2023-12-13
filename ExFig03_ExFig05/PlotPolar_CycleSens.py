#! /usr/bin/env python 


import argparse
import numpy as np
import numpy.random
import copy
import matplotlib.pyplot as plt
import scipy.stats as st
from scipy.stats import bootstrap
from scipy import integrate
from scipy import stats
from astropy import stats as statsAs
from netCDF4 import Dataset
import sys
import matplotlib as mpl
from matplotlib.patches import Rectangle
from tueplots import figsizes


def create_subplot_labels(ax,row,col):
    PLabelMat = [['(a)', '(b)', '(c)', '(d)'],['(e)', '(f)', '(g)', '(h)'],
                 ['(i)', '(j)', '(k)', '(l)',],['(m)', '(n)', '(o)', '(p)'],
                 ['(q)', '(r)', '(s)', '(t)'], ['(u)', '(v)', '(w)', '(x)'],
                 ['(y)', '(z)', '(aa)', '(bb)']]
    ax.text(0.05,0.95, PLabelMat[row][col],color='black',
            horizontalalignment='center', verticalalignment='center',
            transform=ax.transAxes,fontsize=14)

def compute_flux(Time,Flux,Threshold=2.25e12,Clen=1500):
    CLenInd=int(Clen/100)
    # print(Threshold)
    CumFlux=np.zeros(CLenInd)
    CumFluxMax=np.zeros(CLenInd)
    counter=0
    for i in range(0,len(Time),CLenInd):
        # print('OuterLoop:',Time[i],Time[i+len(CumFlux)])
        if (i > CLenInd) and (i+CLenInd < len(Time)) and counter == 0:
            FluxWind=Flux[i:i+len(CumFlux)]
            Con1=(FluxWind>Threshold)*1.0
            MaxVal=FluxWind.max()
            # print('Inner loop:',Time[i],Time[i+len(CumFlux)])
            # sys.exit()
            if MaxVal > Threshold:
                # print(i,i+len(CumFlux))
                # print(Time[i]-65000,Time[i+len(CumFlux)]-65000)
                # print(FluxWind)
                # sys.exit()
                counter=1
                Con2=(FluxWind==MaxVal)*1.0
            else:
                Con2=np.zeros(CLenInd)
            CumFlux=CumFlux+ (FluxWind>Threshold)*1.0
            CumFluxMax=CumFluxMax + Con2
        else:
            counter = 0
    return CumFlux,CumFluxMax

def compute_hist(Time,Flux,Threshold=2.25e12,CLen=1500):
    _,FMax = compute_flux(Time,Flux,Threshold=Threshold,Clen=CLen)
    # print(FMax)
    Event = np.asarray(flux2hetiming(FMax))
    return Event

def get_he_events(Run,Threshold=2.25e12,CLen=1500):
    File=BasePath + Run + "/Postprocessing/" + Region + "Hov.nc"
    T, F = ReadData(File,'flux_crossSection')
    Event = compute_hist(T,F,CLen=CLen, Threshold=Threshold)
    return Event

def extract_surges(icevol,time=100,thresh=-5e11):
    vol_grad = np.gradient(icevol,time)
    # fig,ax1 = plt.subplots(1,1)
    # print(vol_grad)
    # ax1.plot(T,icevol)
    surge_inds = [i for i, n in enumerate(vol_grad) if n < thresh]
    surge_list = [[]] # array of sub-arrays
    for i, num in enumerate(surge_inds):          # go through each element after the first
        surge_list[-1].append(num)                  # Add it to the last sub-array
        if i != len(surge_inds)-1 and (surge_inds[i+1] - surge_inds[i]) > 1 :   # if 5 encountered and not last element
            surge_list.append([])
    filtered_surge_list=[x for x in surge_list if len(x)>=4]
    # print(filtered_surge_list)
    # for j in range(len(filtered_surge_list)):
        # ax1.plot(T[filtered_surge_list[j]],icevol[filtered_surge_list[j]], color='green', lw=3)
    # plt.show()
    # sys.exit()
    return filtered_surge_list

def flux2hetiming(FluxVec,spacing=100):
    HETiming = []
    # print(FluxVec)
    for i in range(len(FluxVec)):
        for j in range(int(FluxVec[i])):
            HETiming.append(i*spacing)
    return HETiming

def kuipers_test(Event,BaseSample,FS=8):
    KTest = abs(round(statsAs.kuiper_two(BaseSample,Event)[1],3))
    KTestStr = "p-value: " + str(KTest)
    return (1-KTest)*100

def plot_do_cycle_polar(ax,CLen=1500):
    n = 200
    if CLen == 1000:
        thetaC1 = np.linspace(0, 2*np.pi / CLen*100, n)
        thetaC2 = np.linspace(2*np.pi / CLen*900,2*np.pi/CLen*CLen, n)
        thetaW1 = np.linspace(2*np.pi / CLen*100,2*np.pi/CLen*500, n)
        thetaW2 = np.linspace(2*np.pi / CLen*500,2*np.pi/CLen*900, n)
    elif CLen == 2000:
        thetaC1 = np.linspace(0, 2*np.pi / CLen*300, n)
        thetaC2 = np.linspace(2*np.pi / CLen*1700,2*np.pi/CLen*CLen, n)
        thetaW1 = np.linspace(2*np.pi / CLen*300,2*np.pi/CLen*1000, n)
        thetaW2 = np.linspace(2*np.pi / CLen*1000,2*np.pi/CLen*1700, n)
    t = np.linspace(0,2*np.pi,n)   #theta values
    r = np.linspace(2.10e15,5.9e15,2)        #radius values full circle
    rg, tg = np.meshgrid(r,thetaC1)      #create a r,theta meshgrid
    c = tg                         #define color values as theta value
    norm = mpl.colors.Normalize(0-1, thetaC1[-1])
    im = ax.pcolormesh(thetaC1, r, c.T, norm=norm, cmap='Reds_r')
    rg, tg = np.meshgrid(r,thetaW1)
    c = tg
    norm = mpl.colors.Normalize(thetaW1[0], thetaW1[-1]+1)
    im = ax.pcolormesh(thetaW1, r, c.T,norm=norm, cmap='Blues')
    rg, tg = np.meshgrid(r,thetaW2)      #create a r,theta meshgrid
    c = tg                         #define color values as theta value
    norm = mpl.colors.Normalize(thetaW2[0]-1, thetaW2[-1])
    im = ax.pcolormesh(thetaW2, r, c.T,norm=norm, cmap='Blues_r')
    rg, tg = np.meshgrid(r,thetaC2)      #create a r,theta meshgrid
    c = tg                         #define color values as theta value
    norm = mpl.colors.Normalize(thetaC2[0], thetaC2[-1]+1)
    im = ax.pcolormesh(thetaC2, r, c.T,norm=norm, cmap='Reds')

def plot_he_lines(ax, Theta, Vol, surge_events):
    ax.plot(Theta,Vol,color="black",lw=0.1)
    for i in range(len(surge_events)):
        if i == 0:
            ax.plot(Theta[surge_events[i]],Vol[surge_events[i]], color='black', lw=1,
                label='Ice volume')
        else:
            ax.plot(Theta[surge_events[i]],Vol[surge_events[i]], color='black', lw=1)

def plot_panel_titles(ax,iSMB,iTemp):
    if iSMB == 0:
        ax[iSMB][iTemp].set_title("$\Delta$Temp = " + str(TPert[iTemp])+"°C", \
            fontsize=14,fontweight='bold')
    if iTemp == 0:
        SMBStr = "$\Delta$SMB = " + str(round(SMBPert[iSMB]/1000,3))+" m/yr"
        ax[iSMB][iTemp].text(-0.2, 0.5, SMBStr,color="black",
        horizontalalignment='center', verticalalignment='center',
        transform=ax[iSMB][iTemp].transAxes, rotation=90,
        fontsize=14,fontweight='bold')

def plot_siglevel(SigLevel,ax):
    Alpha=0.2
    rect  = [0.1,0.1,0.8,0.8]
    axx = ax.inset_axes([.0, .0, 1, 1],)
    autoAxis = axx.axis('off')
    if SigLevel >= 99.0:
        rec = Rectangle((-0.1,-0.05),1.2,1.1,fill=True,
                facecolor='lime',alpha=Alpha)
        rec = axx.add_patch(rec)
        rec.set_clip_on(False)
    elif SigLevel >= 90.0 and SigLevel < 99.0:
        rec = Rectangle((-0.1,-0.05),1.2,1.1,fill=True,
                facecolor='grey',alpha=Alpha)
        rec = axx.add_patch(rec)
        rec.set_clip_on(False)
    axx.set_zorder(0)
    ax.set_zorder(axx.get_zorder()+1)

def ReadData(FilePath,var='thk'):
    NCFile=Dataset(FilePath)
    Time=NCFile.variables['time'][:]/(86400*365)
    Time = Time - Time[0]
    if var == 'thk':
        Vol=NCFile.variables[var][:,0,0]
        return Time,Vol
    else:
        Flux=NCFile.variables['flux_crossSection'][:]
        return Time,Flux
    # [print(FilePath,Time[x]) for x in range(len(Time))]

def set_axis_properties(ax,xticks,labels):
    if Region == "Kenzie":
        ax.set_ylim(3.85e15,4.95e15)
    else:
        ax.set_ylim(4e15,5.9e15)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1) # clockwise
    ax.set_yticklabels([])
    ax.yaxis.grid(False)
    ax.set_xticks(xticks)
    ax.set_xticklabels(labels)
    ax.xaxis.grid(True)
    ax.tick_params(pad=2)


# Ticks=np.array([0, 300, 600, 900, 1200])
# xticks = 2*np.pi/1500.0*Ticks
# labels = ['0', '300', '600', '900', '1200']
# counter = 0

# # fig,ax = plt.subplots(7,4,figsize=(8,14),constrained_layout=True, subplot_kw={'projection': 'polar'})
# m,n = 6,4
# # fig,ax = plt.subplots(m,n,subplot_kw={'projection': \
    # # 'polar'},gridspec_kw={'wspace':0.1,'hspace':0.175,'top':0.99, 'bottom':0.01, \
        # # 'left':0.03, 'right':0.97}, figsize=(17,14))
# fig,ax = plt.subplots(m,n,subplot_kw={'projection': \
    # 'polar'},gridspec_kw={'wspace':0.15,'hspace':0.175,'top':0.99, 'bottom':0.01, \
        # 'left':0.03, 'right':0.97}, figsize=(10.5,14))
# Run=Runs[1][2]

def Main():
    global BasePath, Region, TPert, SMBPert
    parser = argparse.ArgumentParser()
    parser.add_argument("-reg","--region", choices=["Hudson",
        "Kenzie"], help="Specify which region you want to plot")
    parser.add_argument("-clen","--cyclelen", choices=["1000",
        "2000"], help="Specify DO cycle length")

    args = parser.parse_args()
    Region = args.region
    CLen = float(args.cyclelen)
    TPert = [14,10,6,2]
    SMBPert = [100,50]
    BasePath="/work/ba0989/m300792/HE_Runs/ExperimentsComposite/"
    if args.region == "Kenzie":
        BaseRunSample = get_he_events("HE02",Threshold=1.0e12, CLen=CLen)
        GradThresh = -2e11
        FluxThresh = 1.5e12
    else:
        BaseRunSample = get_he_events("HE02", CLen=CLen)
        GradThresh = -3e11
        FluxThresh = 2.25e12
    if args.cyclelen == "1000":
        Runs = [["HE81", "HE86"], ["HE78", "HE87"], ["HE79", "HE88"],\
                   ["HE80", "HE89"]]
        Ticks=np.array([0, 250, 500, 750])
        labels = ['0', '250', '500', '750']
    elif args.cyclelen == "2000":
        Runs = [["HE77", "HE82"], ["HE74", "HE83"], ["HE75", "HE84"],\
                   ["HE76", "HE85"]]
        Ticks=np.array([0, 500, 1000, 1500])
        labels = ['0', '500', '1000', '1500']
        print(BaseRunSample)
        # sys.exit()
    xticks = 2*np.pi/CLen*Ticks
    counter = 0
    m,n = 4,2
    fig,ax = plt.subplots(n,m,subplot_kw={'projection': \
        'polar'},gridspec_kw={'wspace':0.43,'hspace':0.175,'top':0.99, 'bottom':0.01, \
            'left':0.03, 'right':0.97}, figsize=(12,5))

    for iSMB in range(len(Runs[0])):
        for iTemp in range(len(Runs)):
            File=BasePath + Runs[iTemp][iSMB] + "/Postprocessing/IceVolume" +\
            Region + ".nc"
            global T
            T, Vol = ReadData(File)
            HESample = get_he_events(Runs[iTemp][iSMB], Threshold=FluxThresh, CLen=CLen)
            sig_level = kuipers_test(HESample,BaseRunSample)
            plot_siglevel(sig_level,ax[iSMB][iTemp])
            print(iTemp, iSMB,Runs[iTemp][iSMB], sig_level)
            print(HESample)
            surge_events = extract_surges(Vol,thresh=GradThresh)
            Theta = 2*np.pi/CLen*T
            plot_he_lines(ax[iSMB][iTemp], Theta, Vol, surge_events)
            plot_do_cycle_polar(ax[iSMB][iTemp],CLen=CLen)
            set_axis_properties(ax[iSMB][iTemp],xticks,labels)
            plot_panel_titles(ax,iSMB,iTemp)
            create_subplot_labels(ax[iSMB][iTemp],iSMB,iTemp)

    ax[0,0].legend(loc=1, bbox_to_anchor=(1.25,1.1))
    plt.savefig('Polar_' + args.cyclelen + '_' + Region + '.png', dpi=300, bbox_inches='tight')
    plt.show()
    sys.exit()

if __name__ == '__main__':
    Main()
