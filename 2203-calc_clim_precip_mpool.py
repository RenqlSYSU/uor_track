#!/usr/bin/env python
'''
1. when cyclone located in the defined region, read datetime
2. then convert hourly date to 6 hourly date,
3. composite corresponding u,v,t,z and their variances 
4. store the data

20211014
renql
'''
import xarray as xr
import numpy as np
import pandas as pd
from multiprocessing import Pool
import sys
from datetime import datetime
from renql import cyc_filter
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmaps

def draw_precip(var,cnlev,figtitle,cblabel,figdir):
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
    line = ff.readline()
    while line:
        term = line.strip().split(" ")
        if term[0] == "TRACK_ID":
            linenum = ff.readline()
            term1 =linenum.strip().split(" ")
            num = int(term1[-1])
            
            for nl in range(0,num,1):
                line = ff.readline()
                data = list(map(float,line.strip().split(" ")))
                if data[1]<=flonr and data[1] >= flonl and\
                data[2]<=flatn and data[2]>=flats :
                    ctime.append(datetime.strptime(str(int(data[0])),'%Y%m%d%H'))
                    clat.append(data[2])
                    clon.append(data[1])

        line = ff.readline()
    ff.close()
    ctime = xr.DataArray(ctime)
    clat = xr.DataArray(clat)
    clon = xr.DataArray(clon)
    return ctime,clat,clon

def mask_precip(nl):
    ds  = xr.open_dataset(datapath+"1980.nc")
    lat = ds.latitude
    lon = ds.longitude
    ilon = lon[(lon>=lonl) & (lon<=lonr)]
    ilat = lat[(lat>=lats) & (lat<=latn)]
    '''
    term = ds['tp'].sel(time=ds.time.dt.year.isin(1980)
            ,longitude=ilon,latitude=ilat)
    var = term.groupby(term.time.dt.month).sum('time') 
    for ny in range(1981,2021,1):
        print(ny)
        ds  = xr.open_dataset("%s%d.nc"%(datapath,ny))
        term = ds['tp'].sel(time=ds.time.dt.year.isin(ny)
                ,longitude=ilon,latitude=ilat)
        var = var + term.groupby(term.time.dt.month).sum('time') 
    var = var*1000/41
    var.attrs['units']='mm/month'
    ds1 = var.to_dataset(name='tp')
    ds1.to_netcdf("%sclim_precip.nc"%(fileout),"w")
    '''
    ctime,clat,clon = composite_time("%s%s_%d_1980-2020_%s"%(
        path,prefix,nl,suffix),lats,latn,lonl,lonr)

    var  = np.empty( [12,len(ilat),len(ilon)],dtype=float )  
    for ny in range(1980,2021,1):
        print('task [%d]:%d'%(nl,ny))
        ds  = xr.open_dataset("%s%d.nc"%(datapath,ny))
        term = ds['tp'].sel(time=ds.time.dt.year.isin(ny)
                ,longitude=ilon,latitude=ilat)
        
        ctime1 = ctime.where(ctime.dt.year.isin(ny),drop=True)
        clat1 = clat.where(ctime.dt.year.isin(ny),drop=True)
        clon1 = clon.where(ctime.dt.year.isin(ny),drop=True)
        for ct,clo,cla in zip(ctime1,clon1,clat1):
            indx = np.argwhere(term.time.data==ct.data)[0][0]
            term[indx:(indx+3),:,:] = term[indx:(indx+3),:,:].where(
                (np.square(term.longitude-clo)+np.square(term.latitude-cla))>(radiu*radiu), 0)
        #term = term/term.time.dt.days_in_month
        var = var + term.groupby(term.time.dt.month).sum('time') 
    var = var*1000/41
    var.attrs['units']='mm/month'
    ds1 = var.to_dataset(name='tp')
    ds1.to_netcdf("%sclim_precip_%d_%s.nc"%(fileout,nl,suffix),"w")

def mask_precip_oneyear(ny):
    ds  = xr.open_dataset("%s%d.nc"%(datapath,ny))
    lat = ds.latitude
    lon = ds.longitude
    ilon = lon[(lon>=lonl) & (lon<=lonr)]
    ilat = lat[(lat>=lats) & (lat<=latn)]
    term = ds['tp'].sel(time=ds.time.dt.year.isin(ny)
            ,longitude=ilon,latitude=ilat)
    var = term.groupby(term.time.dt.month).sum('time') 
    var.attrs['units']='mm/month'
    ds1 = var.to_dataset(name='tp')
    ds1.to_netcdf("%s%d_precip.nc"%(fileout,ny),"w")

    for nl in lev:
        ctime,clat,clon = composite_time("%s%s_%d_1980-2020_%s"%(
            path,prefix,nl,suffix),lats,latn,lonl,lonr)

        ctime1 = ctime.where(ctime.dt.year.isin(ny),drop=True)
        clat1 = clat.where(ctime.dt.year.isin(ny),drop=True)
        clon1 = clon.where(ctime.dt.year.isin(ny),drop=True)
        print(ctime1)
        for ct,clo,cla in zip(ctime1,clon1,clat1):
            print(ct)
            term.loc[ct,:,:] = term.sel(time=ct).where(
                (np.square(term.longitude-clo)+np.square(term.latitude-cla))>25, 0)
        var.data = term.groupby(term.time.dt.month).sum('time').data
        ds1 = var.to_dataset(name='tp')
        ds1.to_netcdf("%s%d_precip_%d_%s.nc"%(fileout,ny,nl,suffix),"w")

flats = 25  #int(sys.argv[2])
flatn = 45  #int(sys.argv[3])
flonl = 60  #int(sys.argv[4])
flonr = 105 #int(sys.argv[5])
prefix = "fft"
#suffix = '0_2545-60110'
#suffix = '5_2545-60110_2_4545-60110'
suffix = '5_2545-60110_2_2545-6060'
behv = ["ALL" ,"PAS" ,"NTP" ,"STP" ,"NTL" ,"STL" ,"LYS" ]#,"DIF"]
title= {'0_2545-60110':'local',
        '5_2545-60110_2_4545-60110':'northern',
        '5_2545-60110_2_2545-6060':'western'}
radiu=5

fileout="/home/users/qd201969/uor_track/mdata/"
lev  = [850,500,250]
path = '/home/users/qd201969/ERA5-1HR-lev/'
figdir = '/home/users/qd201969/uor_track/fig/'
datapath = "/gws/nopw/j04/ncas_generic/users/renql/ERA5_hourly/precip/ERA5_precip_1hr_dec-jan" #1980.nc

lonl=0  #0  #
lonr=150#360#
lats=15  #
latn=70 #
sy='clim'
'''
process_pool = Pool(processes=3)
results=[]
for nl in range(len(lev)):
    result=process_pool.apply_async(mask_precip, args=(lev[nl],))
    results.append(result)
print(results) 
print(results[0].get()) 
process_pool.close()
process_pool.join() 
print(results[0].get()) 
#mask_precip_oneyear(sy)
#mask_precip()
'''
ds = xr.open_dataset("%s%s_precip.nc"%(fileout,sy))
ilon = ds.longitude
ilat = ds.latitude
var = ds['tp']
#draw_precip(var,[0,255,15],'%s'%sy,'precip (mm/month)',figdir+'clim_precip')

for nl in lev:
    ds = xr.open_dataset("%s%s_precip_%d_%s.nc"%(fileout,sy,nl,suffix))
    #ds = xr.open_dataset("%s%s_precip_%d_%s_%ddegree.nc"%(fileout,sy,nl,suffix,radiu))
    term = (var-ds['tp'].data)*100/var
    draw_precip(term,[2,104,6],'%s %s %d'%(title[suffix],sy,nl),
            'precip percent','%sclim_precip%d_%s'%(figdir,nl,suffix))
