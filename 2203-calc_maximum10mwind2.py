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

flats = 25  #int(sys.argv[2])
flatn = 45  #int(sys.argv[3])
flonl = 60  #int(sys.argv[4])
flonr = 105 #int(sys.argv[5])
prefix = "fft"
#suffix = '0_2545-60110'
#suffix = '5_2545-60110_2_4545-60110'
suffix = '5_2545-60110_2_2545-6060'
title= {'0_2545-60110':'local',
        '5_2545-60110_2_4545-60110':'northern',
        '5_2545-60110_2_2545-6060':'western'}
radiu=5
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

def main_run():
    perc = 99
    #max10mwind_threshold(perc)
    ds = xr.open_dataset("%smax10mwind_%dthreshold_month.nc"%(fileout,perc))
    ilon = ds.lon
    ilat = ds.lat
    thre = ds['threshold'].data
    #monthly_contour(thre,ilon,ilat,[2,19,1],'max10mwind','m/s',figdir+'max10mwind')
    
    var  = np.empty( [12,len(ilat),len(ilon)],dtype=float ) 
    #ctrl_max10mwind_event(thre,var)
    #count_max10mwind_event(thre,var,850)
    process_pool = Pool(processes=3)
    results=[]
    for nl in range(len(lev)):
        result=process_pool.apply_async(count_max10mwind_event,
                args=(perc,thre,var,lev[nl],))
        results.append(result)
    print(results) 
    process_pool.close()
    process_pool.join() 
    print(results[0].get()) 
    
    '''
    days =[31   ,28   ,31   ,30   ,31   ,30   ,31   ,31   ,30   ,31   ,30   ,31   ]
    ds = xr.open_dataset("%sclim_max10mwind_event.nc"%(fileout))
    ilon = ds.lon
    ilat = ds.lat
    var = ds['event'].data
    #term = var
    #for i in range(len(days)):
    #    term[i,:,:] = var[i,:,:]/days[i]
    #monthly_contour(term*30,ilon,ilat,[23,40,1],'max10mwind','h/30day',figdir+'max10mwind')

    for nl in lev:
        ds = xr.open_dataset("%sclim_max10mwind_event_%d_%s.nc"%(fileout,nl,suffix))
        term = (var-ds['event'].data)*100/var
        monthly_contour(term,ilon,ilat,[2,104,6],'%s %d'%(
            title[suffix],nl),'max10mwind (%)','%smax10mwind%d_%s'%(figdir,nl,suffix))
    '''
def max10mwind_threshold(perc):
    ds  = xr.open_dataset(datapath+"1980.nc")
    ilon = ds.lon
    ilat = ds.lat
    thre = np.empty( [12,len(ilat),len(ilon)],dtype=float ) 
    
    for nm in range(0,12,1):
        var = ds['var1'].sel(time=np.array(
            [ds.time.dt.month.isin(nm+1),ds.time.dt.year.isin(1980)]
            ).all(axis=0)).data
        for ny in range(1981,2021,1):
            ds1  = xr.open_dataset("%s%d.nc"%(datapath,ny))
            term = ds1['var1'].sel(time=np.array(
                [ds1.time.dt.month.isin(nm+1),ds1.time.dt.year.isin(ny)]
                ).all(axis=0)).data
            var = np.concatenate((var, term))
            print("month %2d %d: "%(nm,ny), var.shape)
        thre[nm,:,:] = np.percentile(var,perc,axis=0)
    
    da = xr.DataArray(thre, coords=[range(0,12,1),ilat,ilon], 
            dims=["month","lat","lon"])
    ds2 = da.to_dataset(name='threshold')
    ds2.to_netcdf("%smax10mwind_%dthreshold_month.nc"%(fileout,perc),"w")

def count_max10mwind_event(perc,thre,var,lev):
    ctime,clat,clon = composite_time("%s%s_%d_1980-2020_%s"%(
        path,prefix,lev,suffix),lats,latn,lonl,lonr)

    for ny in range(1980,2021,1):
        print('task [%d]:%d'%(lev,ny))
        ds  = xr.open_dataset("%s%d.nc"%(datapath,ny))
        term = ds['var1'].sel(time=ds.time.dt.year.isin(ny))
        for nm in range(12):
            term.loc[dict(time=term.time.dt.month.isin(nm+1))] = xr.where(
                term.sel(time=term.time.dt.month.isin(nm+1))>thre[nm,:,:],1,0)
        
        ctime1 = ctime.where(ctime.dt.year.isin(ny),drop=True)
        clat1 = clat.where(ctime.dt.year.isin(ny),drop=True)
        clon1 = clon.where(ctime.dt.year.isin(ny),drop=True)
        for ct,clo,cla in zip(ctime1,clon1,clat1):
            if ct.dt.dayofyear==366:
                continue
            #indx = np.argwhere(term.time.data==ct.data)[0][0]
            #term[indx,:,:] = term[indx,:,:].where(
            term.loc[ct,:,:] = term.sel(time=ct).where(
                (np.square(term.lon-clo)+np.square(term.lat-cla))>(radiu*radiu), 0)
        var = var + term.groupby(term.time.dt.month).sum('time') 
    var = var/41
    ds1 = var.to_dataset(name='event')
    ds1.to_netcdf("%sclim_%dmax10mwind_event_%d_%s.nc"%(fileout,perc,lev,suffix),"w")

def composite_time(filname,flats,flatn,flonl,flonr):
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
                #data = list(map(float,line.strip().replace(" &","").split(" ")))
                #if data[0]<=9504 and data[0]>=745 and \
                #data[7]<=flonr and data[7]>=flonl and\
                #data[8]<=flatn and data[8]>=flats :
                data = list(map(float,line.strip().split(" ")))
                if data[1]<=flonr and data[1] >= flonl and\
                data[2]<=flatn and data[2]>=flats :
                    ctime.append(datetime.strptime(str(int(data[0])),'%Y%m%d%H'))
                    clat.append(data[2])
                    clon.append(data[1])
                    #clon.append(data[7])
                    #clat.append(data[8])
                    #mswd.append(data[9])

        line = ff.readline()
    ff.close()
    ctime = xr.DataArray(ctime)
    clat = xr.DataArray(clat)
    clon = xr.DataArray(clon)
    return ctime,clat,clon

def monthly_contour(var,ilon,ilat,cnlev,figtitle,cblabel,figdir):
    print(var.min())
    print(var.max())
    lat_sp = 20
    lon_sp = 30 #60 #
    nrow = 4 #6 #
    ncol = 3 #2 #
    bmlo = 0.37 #0.25 #
    title_font=14
    label_font=10
    titls=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    
    ds = xr.open_dataset('/gws/nopw/j04/ncas_generic/users/renql/ERA5_mon/ERA5_mon_u_1979-2020.nc')
    da = ds['u'].sel(level=200,longitude=ilon,latitude=ilat,method="nearest").load()
    uwnd = da.groupby(da.time.dt.month).mean('time')
    del ds, da

    ds = xr.open_dataset("/home/users/qd201969/gtopo30_0.9x1.25.nc")
    phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").load()
    phis = phis/9.8 # transfer from m2/s2 to m
    del ds
    
    cnlevels = np.arange(cnlev[0], cnlev[1], cnlev[2])
    fcolors = cmaps.precip2_17lev
    norm = colors.BoundaryNorm(boundaries=cnlevels, 
        ncolors=fcolors.N,extend='both')
    
    fig = plt.figure(figsize=(12,12),dpi=300)
    ax = fig.subplots(nrow, ncol, subplot_kw=dict(
        projection=ccrs.PlateCarree(central_longitude=180.0)))
    nm = -1
    for nr in range(0,nrow,1):
        for nc in range(0,ncol,1):
            nm = nm+1 
            axe = ax[nr][nc]
            axe.add_feature(cfeat.GSHHSFeature(levels=[1,2],edgecolor='k')
                    , linewidth=0.8, zorder=1)
            axe.set_title(figtitle+" "+titls[nm],fontsize=title_font)

            cont = axe.contourf(ilon, ilat, var[nm,:,:], cnlevels, 
                 transform=ccrs.PlateCarree(),cmap=fcolors,
                 extend='both',norm=norm)
            topo = axe.contour(ilon, ilat, phis, [1500,3000], 
                 transform=ccrs.PlateCarree(),colors='black',linewidths=1.5)
            axe.plot([flonl,flonl,flonr,flonr,flonl],[flatn,flats,flats,flatn,flatn], 
                 linewidth=2.5, color='black', transform=ccrs.PlateCarree()) # filter box
            jets = axe.contour(ilon, ilat, uwnd[nm,:,:], [30,40,50], 
                 transform=ccrs.PlateCarree(),colors='darkviolet',linewidths=1.5)
            
            if nc == 0:
                axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
                axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
            if nr == (nrow-1):
                axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
                axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))

    if (nrow/ncol) >= 2.0: 
        position = fig.add_axes([0.99, bmlo+0.05, 0.01, 0.6]) #left, bottom, width, height
        cb = plt.colorbar(cont, cax=position ,orientation='vertical')#, shrink=.9)
        cb.set_label(label=cblabel, size=title_font) #, weight='bold'
    else:
        position = fig.add_axes([0.2, bmlo+0.005, 0.7, 0.01]) #left, bottom, width, height
        cb = plt.colorbar(cont, cax=position ,orientation='horizontal')#, shrink=.9)
        plt.figtext(0.02,bmlo-0.005, cblabel,fontsize=title_font,
                horizontalalignment='left',verticalalignment='bottom')
    
    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig(figdir, bbox_inches='tight',pad_inches=0.01)

if __name__=='__main__':
    main_run()
