#!/usr/bin/env python

import sys, os, subprocess, linecache
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmaps

'''
read total cyclone in ff_250_1980-2020_2_3045-5960
get the lifetime, max intensity

20210928
'''
def behavior(filname,flats,flatn,flonl,flonr):
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

    #===============================================
    # draw figure 
    #===================================================
    lonl=0  #0  #
    lonr=150#360#
    lats=15 #0  #
    latn=70 #90 #
    lat_sp = 20
    lon_sp = 30
    nrow = 2
    ncol = 1
    bmlo = 0.4
    BIGFONT=22
    MIDFONT=14
    SMFONT=10

    ds = xr.open_dataset("/home/users/qd201969/gtopo30_0.9x1.25.nc")
    lat = ds.lat
    lon = ds.lon
    ilon = lon[(lon>=lonl) & (lon<=lonr)]
    ilat = lat[(lat>=lats) & (lat<=latn)]
    phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").load()
    phis = phis/9.8 # transfer from m2/s2 to m
    del ds

    term = filname.split("/")
    fig = plt.figure(figsize=(9,9),dpi=200)

    axe = plt.subplot(nrow,ncol,1,projection=ccrs.PlateCarree())    #创建子图
    axe.add_feature(cfeat.GSHHSFeature(levels=[1,2],edgecolor='k'), linewidth=0.8, zorder=1)
    axe.set_title(term[-1],fontsize=MIDFONT)

    axe.contour(ilon, ilat, phis, [1500,3000], 
          transform=ccrs.PlateCarree(),colors='black',linewidths=1.5)
    axe.scatter(tlon,tlat,transform=ccrs.PlateCarree())
    axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
    axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
    axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
    axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))

    ax = plt.subplot(nrow,ncol,2)    #创建子图
    ax.scatter(life,inte)
    ax.set_xlabel("lifetime (days)",fontsize=MIDFONT)
    ax.set_ylabel("Max intensity ($10^{-5} s^{-1}$)",fontsize=MIDFONT)
    ax.set_title(term[-1],fontsize=MIDFONT)

    fig.savefig("/home/users/qd201969/uor_track/fig/life_inte_"+term[-1]+".png")

def monthly_contour(var,cnlev,figtitle,cblabel,figdir):
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

