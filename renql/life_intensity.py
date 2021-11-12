#!/usr/bin/env python
'''
read total cyclone number in ff_250_1980-2020_2_3045-5960
then use box to filter different behaviors cyclones
plot the filter box in the trajectory figure

20210928
'''

import sys, os, subprocess, linecache
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib import colors
import cmaps

def calc_life_intensity(filname):
    ff = open(filname,"r") 
    line1 = ff.readline()
    line2 = ff.readline()
    line3 = ff.readline()
    line4 = ff.readline()

    tid = []
    life = []
    inte = [] #  max intensity
    tlat = []  # max intensity location
    tlon = []  # max intensity location
    line = ff.readline()
    while line:
        term = line.strip().split(" ")
        if term[0] == "TRACK_ID":
            tid.append(term[-1])
            
            linenum = ff.readline()
            term1 =linenum.strip().split(" ")
            num = int(term1[-1])
            life.append(num/24.0)
            
            data=[]
            for nl in range(0,num,1):
                data.append(list(map(float, ff.readline().strip().split(" "))))

            data = np.array(data)
            inte.append(data[:,3].max())
            loc = np.argmax(data[:,3])
            tlat.append(data[loc,2])
            tlon.append(data[loc,1])

        line = ff.readline()

    a = line4.strip().split(" ",1)
    term = a[1].strip().split(" ",1)
    print("total cyclone number in %s : %s" %(ff.name,term[0]))
    ff.close()

    print("tid life(days) intensity longitude latitude")
    for ni in range(0,len(tid),1):
        print("%s "%tid[ni] + "%f "*4 %(life[ni],inte[ni],tlon[ni],tlat[ni]))

    return life, inte

def hist_life_intensity(filname):
    BIGFONT=22
    MIDFONT=14
    SMFONT=10

    cnlevels = np.arange(0,85,5) #500Z, gpm
    norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=cmaps.precip2_17lev.N,extend='both')
    term = filname.split("/")
    xbin = range(1,16,1)
    ybin = range(1,21,1)
    
    fig = plt.figure(figsize=(9,9),dpi=200)
    axe = fig.subplots(nrow, ncol, sharex=True, sharey=True)

    ax = axe[nr][nc]
    h,xedge,yedge,patches = ax.hist2d(life,inte,bins=[xbin,ybin],weight=True,norm=norm)
    print(h,xedge,yedge)
    ax.set_xlabel("lifetime (days)",fontsize=MIDFONT)
    ax.set_ylabel("Max intensity ($10^{-5} s^{-1}$)",fontsize=MIDFONT)
    ax.set_title(term[-1],fontsize=MIDFONT)

    plt.colorbar(patches, ax=ax, shrink=.7)
    fig.savefig("/home/users/qd201969/uor_track/fig/life_inte_"+term[-1]+".png")

