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

def calc_life_intensity(filname,flats=25,flatn=48,flonl=57,flonr=100):
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
                term = list(map(float, ff.readline().strip().split(" ")))
                if term[1]>=flonl and term[1]<=flonr and term[2]>=flats and term[2]<=flatn:
                    data.append(term)

            data = np.array(data)
            if len(data.shape) == 2:
                inte.append(data[:,3].mean())
                #inte.append(data[:,3].max())
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

    return life, inte, int(term[0])

def hist_life_intensity(filname,ax=None,title=None,figsave=False,flats=25,flatn=45,flonl=60,flonr=90):
    title_font=18
    label_font=14

    life, inte, numb = calc_life_intensity(filname,\
        flats=flats-3,flatn=flatn+3,flonl=flonl-3,flonr=flonr+3)

    #cnlevels = np.arange(0,85,5)
    cnlevels = np.arange(0.5,17.5,1)/100.0
    norm  = colors.BoundaryNorm(boundaries=cnlevels, ncolors=cmaps.precip2_17lev.N,extend='both')
    xbin = np.arange(1,8,0.5)
    ybin = np.arange(1,13,0.5)
   
    if not title:
        title = filname.split("/")[-1]
    
    if not ax:
        fig = plt.figure(figsize=(9,9),dpi=200)
        ax = plt.axes()
        ax.set_xlabel("lifetime (days)",fontsize=label_font)
        ax.set_ylabel("Max intensity ($10^{-5} s^{-1}$)",fontsize=label_font)

    h,xedge,yedge,patches = ax.hist2d(life,inte,bins=[xbin,ybin],\
            cmap=cmaps.precip2_17lev,norm=norm, density=True )
    ax.set_title("%s (%d)"%(title,numb),fontsize=title_font)

    if figsave == True:
        plt.colorbar(patches, ax=ax, shrink=.7)
        fig.savefig("/home/users/qd201969/uor_track/fig/life_inte_"+term[-1]+".png")

    return patches

