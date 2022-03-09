#!/usr/bin/env python
import xarray as xr
import numpy as np
import pandas as pd
from multiprocessing import Pool
import sys, os, subprocess, linecache, gc
from datetime import datetime
from renql import cyc_filter, draw
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmaps

def composite_time(filname,flats,flatn,flonl,flonr):
# prefix of filname must be ffadd
# read max 10m wind as well as its location and time
    ff = open(filname,"r") 
    line1 = ff.readline()
    line2 = ff.readline()
    line3 = ff.readline()
    line4 = ff.readline()
    a = line4.strip().split(" ",1)
    term = a[1].strip().split(" ",1)
    print("total cyclone number in %s : %s" %(ff.name,term[0]))
    
    ctime=[]
    clat =[]
    clon =[]
    mswd =[]
    line = ff.readline()
    while line:
        term = line.strip().split(" ")
        if term[0] == "TRACK_ID":
            linenum = ff.readline()
            term1 =linenum.strip().split(" ")
            num = int(term1[-1])
            
            for nl in range(0,num,1):
                line = ff.readline()
                data = list(map(float,line.strip().replace(" &","").split(" ")))
                if data[0]<=9504 and data[0]>=745 and \
                data[7]<=flonr and data[7]>=flonl and\
                data[8]<=flatn and data[8]>=flats :
                    ctime.append(int(data[0]))
                    clon.append(data[7])
                    clat.append(data[8])
                    mswd.append(data[9])

        line = ff.readline()
    ff.close()
    return ctime,clat,clon,mswd

def max10mwind_percent(i,var,ds,ilon,ilat):
    print(var.time)
    for dtime in var.time:
        indx = np.argwhere(ds.time.data==dtime.data)[0][0]
        indx2=np.arange(indx-24*lagd,indx+24*(lagd+1),24)
        indx3=np.array([np.arange(d-lagh,d+lagh+1) 
            for d in indx2]).reshape((2*lagh+1)*(2*lagd+1))
        
        ds1  = xr.open_dataset("%s1980.nc"%(datapath))
        term = ds1['var1'][indx3].sel(lon=ilon,lat=ilat).data 
        for ny in range(1981,2021,1):
            ds1  = xr.open_dataset("%s%d.nc"%(datapath,ny))
            term = np.concatenate((term,
                ds1['var1'][indx3].sel(lon=ilon,lat=ilat).data))
        print("task %d: %s"%(i,dtime.data),term.shape)
        var.loc[dtime,:,:] = np.percentile(term,95,axis=0)
    return var.data

def max10mwind_threshold():
    ds  = xr.open_dataset(datapath+"1981.nc")
    lat = ds.lat
    lon = ds.lon
    ilon = lon[(lon>=lonl) & (lon<=lonr)]
    ilat = lat[(lat>=lats) & (lat<=latn)]
    var = ds['var1'].sel(time=ds.time.dt.year.isin(1981)
            ,lon=ilon,lat=ilat).load()
    
    ntask = 20
    intvl = len(var.time)//ntask
    process_pool = Pool(processes=ntask)
    results = [process_pool.apply_async(max10mwind_percent, 
        args=(i,var[intvl*i:intvl*(i+1)-1,:,:],ds,ilon,ilat)) 
        for i in range(ntask)]
    print(results) 
    process_pool.close()
    process_pool.join() 
    print(results[0].get()) 
    for i in range(ntask):
        var[intvl*i:intvl*(i+1)-1,:,:].data = results[i].get()
    
    ds2 = var.to_dataset(name='threshold')
    ds2.to_netcdf("%smax10mwind_threshold_lagd%d_lagh%d.nc"%(fileout,lagd,lagh),"w")

def ctrl_max10mwind():
    for ny in range(1981,2021,1):
        print(ny)
        ds  = xr.open_dataset("%s%d.nc"%(datapath,ny))
        term = ds['tp'].sel(time=ds.time.dt.year.isin(ny)
                ,lon=ilon,lat=ilat)
        var = var + term.groupby(term.time.dt.month).sum('time') 
    ds1 = var.to_dataset(name='')
    ds1.to_netcdf("%sclim_precip.nc"%(fileout),"w")

    for nl in lev:
        ctime,clat,clon = composite_time("%s%s_%d_1980-2020_%s"%(
            path,prefix,nl,suffix),lats,latn,lonl,lonr)

        var  = np.empty( [12,len(ilat),len(ilon)],dtype=float )  
        for ny in range(1980,2021,1):
            print(ny)
            ds  = xr.open_dataset("%s%d.nc"%(datapath,ny))
            term = ds['tp'].sel(time=ds.time.dt.year.isin(ny)
                    ,lon=ilon,lat=ilat)
            
            ctime1 = ctime.where(ctime.dt.year.isin(ny),drop=True)
            clat1 = clat.where(ctime.dt.year.isin(ny),drop=True)
            clon1 = clon.where(ctime.dt.year.isin(ny),drop=True)
            print(ctime1)
            for ct,clo,cla in zip(ctime1,clon1,clat1):
                print(ct)
                term.loc[ct,:,:] = term.sel(time=ct).where(
                    (np.square(term.lon-clo)+np.square(term.lat-cla))>25, 0)
            #term = term/term.time.dt.days_in_month
            var = var + term.groupby(term.time.dt.month).sum('time') 
        var = var*1000/41
        var.attrs['units']='mm/month'
        ds1 = var.to_dataset(name='tp')
        ds1.to_netcdf("%sclim_precip_%d_%s.nc"%(fileout,nl,suffix),"w")

flats = 25  #int(sys.argv[2])
flatn = 45  #int(sys.argv[3])
flonl = 60  #int(sys.argv[4])
flonr = 105 #int(sys.argv[5])
prefix = "fft"
suffix = '5_2545-60110_2_2545-6060'
behv = ["ALL" ,"PAS" ,"NTP" ,"STP" ,"NTL" ,"STL" ,"LYS" ]#,"DIF"]
lagd = 2
lagh = 1

fileout="/home/users/qd201969/uor_track/mdata/"
lev  = [850,500,250]
path = '/home/users/qd201969/ERA5-1HR-lev/'
figdir = '/home/users/qd201969/uor_track/fig/'
datapath = "/work/scratch-pw2/renql/ERA5_hourly/wind10/ERA5_speed10_1hr_dec-jan"

lonl=0  #0  #
lonr=150#360#
lats=15  #
latn=70 #
max10mwind_threshold()
#mask_precip_oneyear(1989)
#mask_precip()
'''
sy='1989'
ds = xr.open_dataset("%s%s_precip.nc"%(fileout,sy))
ilon = ds.lon
ilat = ds.lat
var = ds['tp']
#draw.monthly_contour(var*1000,[0,255,15],'%s'%sy,'precip (mm/month)',figdir+'clim_precip')

for nl in lev:
    ds = xr.open_dataset("%s%s_precip_%d_5_2545-60110_2_2545-6060.nc"%(fileout,sy,nl))
    term = (var-ds['tp'])*100/var
    draw.monthly_contour(term,[2,104,6],'%s %d'%(sy,nl),'precip percent','%sclim_precip%d'%(figdir,nl))
'''
