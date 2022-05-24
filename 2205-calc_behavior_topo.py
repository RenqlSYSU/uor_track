#!/usr/bin/env python
'''
1. read /home/users/qd201969/gtopo30_0.9x1.25.nc
    set the phis larger than 1500m as 1
2. if the 6degree cycle of cyclone center has 
    encounter 1, then set this cyclone as the one
    passing through the TP
3. Using location of the first point to define 
    local and remote cyclones

20220518
renql
'''

import xarray as xr
import numpy as np
import pandas as pd
import sys, os, subprocess, linecache, gc
from datetime import datetime
from scipy import stats
from renql import monthly_calc, life_intensity, cyc_filter
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib import colors

title_font=18
label_font=14
plt.rcParams["font.weight"] = "bold"
font = {'family': 'serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black',
        }

if len(sys.argv) < 2 :
    option=0 #int(sys.argv[1]) #Genesis (0)/Lysis (1)/Passing(2)/Passing Time(3)/All Times(4)
    flats = 25  #int(sys.argv[2])
    flatn = 45  #int(sys.argv[3])
    flonl = 60  #int(sys.argv[4])
    flonr = 105 #int(sys.argv[5])
    time = 24 # threshold, hour
    prefix = "fftadd"
    suffix0 = ""
    season = 0 # 0 monthly, 1 seasonal
else:
    option= int(sys.argv[1]) 
    flats = int(sys.argv[2])
    flatn = int(sys.argv[3])
    flonl = int(sys.argv[4])
    flonr = int(sys.argv[5])
    prefix = int(sys.argv[6])
    season = int(sys.argv[7])
    time = int(sys.argv[8])

suffix=str(option)+"_"+str(flats)+str(flatn)+"-"+str(flonl)+str(flonr)
figdir = "/home/users/qd201969/uor_track/fig/"
fileout="/home/users/qd201969/uor_track/mdata/behv4_month_%dh_%s.nc"%(time,suffix)

flonr2 = 110
flatn2 = 40
flats2 = 30
behv = ["ALL" ,"PAS" ,"NTP" ,"STP" ,"NTL" ,"STL" ,"LYS" ]#,"DIF"]
lats = [flats ,flats2,flatn2,flats2,flatn ,flats ,flats ]
latn = [flatn ,flatn2,flatn2,flats2,flatn ,flats ,flatn ]
lonl = [flonl ,flonr2,flonr2,flonr2,flonl ,flonl ,flonl ]
lonr = [flonr ,flonr2,flonr2,flonr2,flonr2,flonr2,flonr2]
lev = [850,500,250]

if season == 0:
    nday=[365,31,28,31,30,31,30,31,31,30,31,30,31]
    month=["ALL","Jan","Feb","Mar","Apr","May","Jun",\
            "Jul","Aug","Sep","Oct","Nov","Dec"]
    frae=744
else:
    nday=[365,90,92,92,91]
    month=["ALL","DJF","MAM","JJA","SON"]
    frae=0

imonth=range(1,len(month),1)
path = '/home/users/qd201969/ERA5-1HR-lev/'
os.chdir("/home/users/qd201969/uor_track/fig")

def main_run():
    #calcbehv()
    #drawannual()
    #draw_table()
    draw_3var_distri()
    #drawbox()

ds = xr.open_dataset("/home/users/qd201969/gtopo30_0.9x1.25.nc")
phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").load()
phis = phis/9.8 # transfer from m2/s2 to m

def draw_3var_distri():
    #varname = ["close1 (%)","close2 (%)","GeopoHeight (m)"]
    varname = ["10mWind (m/s)","MaxRain (mm/d)","RainTime (%)"]
    nrow = len(lev)
    ncol = 3
    bmlo = 0.4
    #xbin = [np.arange(0,101,5), # num=(end-start)/inv
    #        np.arange(0,101,5), 
    #        np.arange(0,8001,400)]
    xbin = [np.arange(2,22.1,1), # num=(end-start)/inv
            np.arange(0,50.1,2.5), 
            np.arange(0,101,5)]
    var = np.empty( [len(varname),len(behv),len(xbin[0])-1],dtype=float )   
    
    fig = plt.figure(figsize=(9,9),dpi=200)
    ax = fig.subplots(nrow, ncol) #sharex=True, sharey=True
    for nl in range(0,len(lev),1):
        #if nl==0 : # for geopotential
        #    xbin[2]=np.arange(1200,1701,25)
        #elif nl==1 :
        #    xbin[2]=np.arange(5000,6001,50)
        #elif nl==2 :
        #    xbin[2]=np.arange(9500,11001,75)

        numb = []
        filname  = path+prefix+"_"+str(lev[nl])+"_1980-2020"
        for nb in range(1,len(behv),1):
            if nb == 0:
                filname1 = filname
            else:
                suffix = str(option[nb])+"_"+str(flats)+str(flatn)+"-"+str(flonl)+str(flonr)
                filname1 = filname+"_"+suffix
            
            life1, inte1, dist1, numb1 = life_intensity.calc_close_cyclone(
                    filname1,flats=20,flatn=50,flonl=105,flonr=130)
            print("%dhPa %s :%d"%(lev[nl],behv[nb],numb1))
            var[0,nb,:],term = np.histogram(life1,xbin[0])
            var[1,nb,:],term = np.histogram(inte1,xbin[1])
            var[2,nb,:],term = np.histogram(dist1,xbin[2])
            var[:,nb,:] = var[:,nb,:]*100/numb1
            numb.append(numb1)
                
        for nv in range(0,3,1):
            ax[nl][nv].grid(True, which="both", color='grey', linestyle='--', linewidth=1)
            ax[nl][nv].yaxis.set_major_formatter(mtick.FormatStrFormatter('%i'))
            for nb in range(1,len(behv),1):
        #        print(var[nv,nb,:])
                ax[nl][nv].plot(xbin[nv][1:],var[nv,nb,:],linewidth=2)
            
            if nl==2 :
                ax[nl][nv].set_xlabel(varname[nv],fontsize=label_font-4,fontdict=font)
            if nv==0 :
                ax[nl][nv].set_ylabel("%dhPa percent"%lev[nl],fontsize=label_font-4,fontdict=font)

    ax[0][2].legend(behv[1:len(behv)], loc='upper right')
    fig.tight_layout(w_pad=0.5,h_pad=1) #,rect=(0,bmlo,1,1)
    fig.savefig("%sclse_pinte_%dh_total"%(figdir,time))

def draw_ts(var,figdir):
    plinsty = ["solid",(0,(4,2)),(0,(1,2))] # change with lev
    pcolor  = ["k","r","b"] # change with option 

    fig = plt.figure(figsize=(9,9),dpi=200)
    axe = plt.axes()
    axe.set_title("%d-%dE, %d-%dN"%(flonl,flonr,flats,flatn),fontsize=title_font,fontdict=font)
    axe.set_xlim(1, 12)
    axe.set_ylim(0, np.max(var))
    axe.set_xlabel("month",fontsize=label_font,fontdict=font)
    axe.set_ylabel("number of tracks",fontsize=label_font,fontdict=font)
    axe.tick_params(axis='both', which='major', labelsize=label_font)
    for nl in range(len(lev)):
        for nop in range(len(option)):
            axe.plot(months,var[nl,nop,:],linewidth=1.2+nl/2.0,color=pcolor[nop],linestyle=plinsty[nl],
                label=str(lev[nl])+" "+behv[nop]) #
    axe.legend(ncol=3,fontsize=label_font)
    fig.savefig(figdir+"annual_ts.png", bbox_inches='tight',pad_inches=0.01)

if __name__=='__main__':
    main_run()

