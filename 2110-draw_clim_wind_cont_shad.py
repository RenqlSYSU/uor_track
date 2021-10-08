#!/usr/bin/env python
'''
read uwnd, vwnd, z to draw monthly wind (vector) and geopotential (shaded)
Loop through the height (850, 500, 250)
maybe later the shaded variable can be changed for t, PV, dtdy

20211007
'''
import sys
import subprocess
import xarray as xr
import numpy as np
import gc #garbage collector
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmaps

lonl=0  #0  #
lonr=150#360#
lats=15 #0  #
latn=70 #90 #
lat_sp = 20
lon_sp = 30
lev = [850,500,250]

nrow = 4
ncol = 3
bmlo = 0.4
BIGFONT=22
MIDFONT=14
SMFONT=10

filename = 'ff_250_1980-2020_2_2545-6080'
files = '/home/users/qd201969/ERA5-1HR-lev/statistic/'+filename+'_stat.nc'
figtitle = '250_2_2545-6080'

titls=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

f = xr.open_dataset('/home/users/qd201969/data/ERA5_mon_z_1979-2020.nc')
lat = f.latitude
lon = f.longitude
ilon = lon[(lon>=lonl) & (lon<=lonr)]
ilat = lat[(lat>=lats) & (lat<=latn)]

ds = xr.open_dataset("/home/users/qd201969/gtopo30_0.9x1.25.nc")
phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").load()
phis = phis/9.8 # transfer from m2/s2 to m

for nl in range(0,len(lev),1):
    da = f['z'].sel(level=lev[nl],longitude=ilon,latitude=ilat,method="nearest").load()
    var = da.groupby(da.time.dt.month).mean('time')

    ds = xr.open_dataset('/home/users/qd201969/data/ERA5_mon_u_1979-2020.nc')
    da = ds['u'].sel(level=lev[nl],longitude=ilon,latitude=ilat,method="nearest").load()
    uwnd = da.groupby(da.time.dt.month).mean('time')

    ds = xr.open_dataset('/home/users/qd201969/data/ERA5_mon_v_1979-2020.nc')
    da = ds['v'].sel(level=lev[nl],longitude=ilon,latitude=ilat,method="nearest").load()
    vwnd = da.groupby(da.time.dt.month).mean('time')
    del ds, da
    gc.collect()

    if cnlev[nl][0] < 0 :
        fcolors = cmaps.BlueDarkRed18
    else:
        fcolors = cmaps.precip2_17lev
    cnlevels = np.arange(cnlev[nl][0], cnlev[nl][1]+cnlev[nl][2], cnlev[nl][2])
    norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=fcolors.N,extend='both')

    fig = plt.figure(figsize=(12,12),dpi=300)
    ax = fig.subplots(nrow, ncol, subplot_kw=dict(projection=ccrs.PlateCarree())) #sharex=True, sharey=True
    nm = -1
    for nr in range(0,nrow,1):
        for nc in range(0,ncol,1):
            nm = nm+1 
            axe = ax[nr][nc]
            #axe = plt.subplot(4,3,nm+1,projection=ccrs.PlateCarree())    #创建子图
            #axe.add_feature(cfeat.COASTLINE.with_scale('110m'),edgecolor='black', linewidth=0.8, zorder=1) 
            #axe.add_feature(cfeat.LAKES.with_scale('110m'),edgecolor='black', linewidth=0.8, zorder=1) 
            axe.add_feature(cfeat.LAND.with_scale('110m'),edgecolor='black', linewidth=0.8, zorder=1) 
            #axe.add_feature(cfeat.GSHHSFeature(levels=[1,2],edgecolor='k'), linewidth=0.8, zorder=1)
            axe.set_title(figtitle+str(lev[nl])+" "+titls[nm],fontsize=SMFONT)

            shad = axe.contourf(ilon, ilat, var[nm,:,:], cnlevels, 
                         transform=ccrs.PlateCarree(),cmap=fcolors,extend='both',norm=norm)
            wind = axe.quiver(ilon[::q_mis], ilat[::q_mis], uwnd[nm,::q_mis,::q_mis],vwnd[nm,::q_mis,::q_mis],
                    zorder=2, pivot='mid',units='inches',scale=30,scale_units='inches',color="dimgray",
                    width=0.02,headwidth=3,headlength=4.5,transform=ccrs.PlateCarree())
            topo = axe.contour(ilon, ilat, phis, [1500,3000], 
                         transform=ccrs.PlateCarree(),colors='black',linewidths=1.5)
            cont = axe.contour(ilon, ilat, var[nm,:,:], [30,32], 
                         transform=ccrs.PlateCarree(),colors='darkviolet',linewidths=2.5)
            if dbox >= 1 :
                axe.plot([flonl,flonl,flonr,flonr,flonl],[flatn,flats,flats,flatn,flatn], 
                         linewidth=2.5, color='black', transform=ccrs.PlateCarree()) # filter box
            #jets = axe.contour(ilon, ilat, uwnd[nm,:,:], [30,32], 
            #             transform=ccrs.PlateCarree(),colors='darkviolet',linewidths=2.5)

            if nc == 0:
                axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
                axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
            if nr == (nrow-1):
                axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
                axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))

    position = fig.add_axes([0.2, bmlo+0.005, 0.7, 0.01]) #left, bottom, width, height
    cb = plt.colorbar(shad, cax=position ,orientation='horizontal')#, shrink=.9)
    axe.quiverkey(wind, 1.05, 0.97, 10, r'$10 m/s$', labelpos='N',coordinates='axes')

    plt.figtext(0.02,bmlo-0.005, var.long_name, fontsize=MIDFONT,
            horizontalalignment='left',verticalalignment='bottom')
    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig(figdir+var.long_name.replace(" ","")+".png", bbox_inches='tight',pad_inches=0.01)

