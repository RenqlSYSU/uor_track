#!/usr/bin/env python
import xarray as xr
import numpy as np
import pandas as pd
from multiprocessing import Pool
import sys, os, subprocess
from datetime import datetime
from renql import cyc_filter
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmaps

suffixs = ['_0_2545-60110','_5_2545-60110_2_4545-60110',
        '_5_2545-60110_2_2545-6060','_2_2545-60110','']
radiu = 6
perc = 99
fileout="/home/users/qd201969/uor_track/mdata/"
figdir = '/home/users/qd201969/uor_track/fig/'
lev  = [850,500,250]
title= {'_0_2545-60110':'local',
        '_5_2545-60110_2_4545-60110':'northern',
        '_5_2545-60110_2_2545-6060':'western',
        '_2_2545-60110':'passing',
        '':'All'}
lonl=0  #0  #
lonr=150#360#
lats=15  #
latn=70 #
flats = 25  #int(sys.argv[2])
flatn = 45  #int(sys.argv[3])
flonl = 60  #int(sys.argv[4])
flonr = 105 #int(sys.argv[5])

def main_run():
    ds = xr.open_dataset("%s/clim_max%dprecip_event.nc"%(fileout,perc))
    ilon = ds.longitude
    ilat = ds.latitude
    var = ds['tp'].data
    draw_annual_contour3x3(var,'%s/clim_max%dprecip_%drad_lag0'%(fileout,perc,radiu),
            [2,104,6],'precip percent',
            '%s/max%dprecip_contribution_3x3.jpg'%(figdir,perc),ilon,ilat)

def calc_associated_weather():
    for suffix in suffixs: 
        #com = "python ~/uor_track/2203-calc_maximum10mwind2.py %s %d %d"\
        #        %(suffix,radiu,perc)
        com = "python ~/uor_track/2203-calc_clim_precip_mpool.py %s %d %d"\
                %(suffix,radiu,perc)
        ret=subprocess.Popen(com,shell=True)
        ret.wait()

def draw_annual_contour3x3(var,filname,cnlev,cblabel,figdir,ilon,ilat):
    lat_sp = 20
    lon_sp = 30 #60 #
    nrow = 3 #6 #
    ncol = 3 #2 #
    bmlo = 0.5 #0.25 #
    title_font=14
    label_font=10
    
    ds = xr.open_dataset('/gws/nopw/j04/ncas_generic/users/renql/ERA5_mon/ERA5_mon_u_1979-2020.nc')
    da = ds['u'].sel(level=200,longitude=ilon,latitude=ilat,method="nearest").load()
    uwnd = da.mean('time').data
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
    for nr in range(0,nrow,1):
        for nc in range(0,ncol,1):
            ds = xr.open_dataset("%s_%d%s.nc"%(filname,lev[nr],suffixs[nc]))
            term = ds['tp'].data
            term = xr.where(var>0,(var-term)*100/var,0)
            term = np.mean(term,axis=0)
            print(term.min())
            print(term.max())
            axe = ax[nr][nc]
            axe.add_feature(cfeat.GSHHSFeature(levels=[1,2],edgecolor='k')
                    , linewidth=0.8, zorder=1)
            axe.set_title('%dhPa %s'%(lev[nr],title[suffixs[nc]]),fontsize=title_font)

            cont = axe.contourf(ilon, ilat, term, cnlevels, 
                 transform=ccrs.PlateCarree(),cmap=fcolors,
                 extend='both',norm=norm)
            topo = axe.contour(ilon, ilat, phis, [1500,3000], 
                 transform=ccrs.PlateCarree(),colors='black',linewidths=1.5)
            axe.plot([flonl,flonl,flonr,flonr,flonl],[flatn,flats,flats,flatn,flatn], 
                 linewidth=2.5, color='black', transform=ccrs.PlateCarree()) # filter box
            jets = axe.contour(ilon, ilat, uwnd, [30,40,50], 
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
