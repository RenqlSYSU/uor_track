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
from datetime import datetime, timedelta

title_font=18
label_font=14
plt.rcParams["font.weight"] = "bold"
font = {'family': 'serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black',
        }

def circle_distance(lat1,lon1,lat2,lon2):
    # Calculate the distance between two points of a sphere
    # The input parameters are scalar, from 1(start) to 2(end)
    radiu = 6378.388 #the radius of earth, km
    lat1 = lat1*np.pi/180.0
    lat2 = lat2*np.pi/180.0

    if lon1>180:
        dlon = (lon2-lon1+360)*np.pi/180.0
    else:
        dlon = (lon2-lon1)*np.pi/180.0
    
    term = np.cos(lat1)*np.cos(lat2)*np.cos(dlon)+\
        np.sin(lat1)*np.sin(lat2)
    dist = radiu*np.arccos(term)

    if abs(dlon) > np.pi:
        dist = 2*radiu*np.pi-dist

    return dist

def calc_life_intensity(filname,flats=25,flatn=45,flonl=57,flonr=110):
    ff = open(filname,"r") 
    line1 = ff.readline()
    line2 = ff.readline()
    line3 = ff.readline()
    line4 = ff.readline()
    
    tid = []
    life = []
    inte = [] #  max intensity
    dist = []
    tlat = []  # max intensity location
    tlon = []  # max intensity location
    line = ff.readline()
    start=datetime(1995, 11, 30, 23, 00)
    while line:
        term = line.strip().split(" ")
        if term[0] == "TRACK_ID":
            tid.append(term[-1])
            
            linenum = ff.readline()
            term1 =linenum.strip().split(" ")
            num = int(term1[-1])
            
            data=[]
            ct1=[]
            for nl in range(0,num,1):
                term = list(map(float, ff.readline().strip().split(" ")))
                ct1.append(start+timedelta(hours=int(term[0])))
                if term[1]>=flonl and term[1]<=flonr and term[2]>=flats and term[2]<=flatn:
                    data.append(term)
                if nl==0:
                    lon1 = term[1]
                    lat1 = term[2]
                if nl==num-1:
                    lon2 = term[1]
                    lat2 = term[2]

            data = np.array(data)
            if sum(i.year==1996 for i in ct1)/len(ct1) >= 0.5 : 
                life.append(num/24.0)
                inte.append(data[:,3].mean())
                #inte.append(data[:,3].max())
                dist.append(circle_distance(lat1,lon1,lat2,lon2))
                loc = np.argmax(data[:,3])
                tlat.append(data[loc,2])
                tlon.append(data[loc,1])

        line = ff.readline()

    a = line4.strip().split(" ",1)
    term = a[1].strip().split(" ",1)
    print("total cyclone number in %s : %s" %(ff.name,term[0]))
    ff.close()

    #print("tid life(days) intensity longitude latitude")
    #for ni in range(0,len(tid),1):
    #    print("%s "%tid[ni] + "%f "*4 %(life[ni],inte[ni],tlon[ni],tlat[ni]))

    return life, inte, dist, len(life)

def hist_life_intensity(filname,ax=None,title=None,figsave=False,\
        flats=25,flatn=45,flonl=60,flonr=90):
    life, inte, dist, numb = calc_life_intensity(filname,\
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
        ax.set_xlabel("lifetime (days)",fontsize=label_font,fontdict=font)
        ax.set_ylabel("Mean intensity ($10^{-5} s^{-1}$)",fontsize=label_font,fontdict=font)

    h,xedge,yedge,patches = ax.hist2d(life,inte,bins=[xbin,ybin],\
            cmap=cmaps.precip2_17lev,norm=norm, density=True )
    ax.set_title("%s (%d)"%(title,numb),fontsize=title_font,fontdict=font)

    if figsave == True:
        plt.colorbar(patches, ax=ax, shrink=.7)
        fig.savefig("/home/users/qd201969/uor_track/fig/life_inte_"+term[-1]+".png")

    return patches

