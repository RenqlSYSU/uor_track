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
from renql import dynamic_calc

lonl=0  #0  #
lonr=150#360#
lats=15 #0  #
latn=70 #90 #
lat_sp = 20
lon_sp = 30
lev = [850,500,250]

title_font=14
label_font=10
plt.rcParams["font.weight"] = "bold"
font = {'family': 'sans-serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black', 
        }
nrow = 4
ncol = 3
bmlo = 0.4

titls= ['DJF','MAM','JJA','SON']
varname = ['vor']
drawvar = ['vor']
unit    = ['s-1']
vcref =[10,20,30] # different levels 
cnlvl =[-3.5,0.5]
cnlvl2=[20,60,100]
q_mis=15
figdir = "/home/users/qd201969/uor_track/fig/"
path = '/gws/nopw/j04/ncas_generic/users/renql/ERA5_mon'

f = xr.open_dataset('%s/ERA5_mon_z_1979-2020.nc'%path)
lat = f.latitude.data
lon = f.longitude.data
ilon = lon[(lon>=lonl) & (lon<=lonr)]
ilat = lat[(lat>=lats) & (lat<=latn)]

ds = xr.open_dataset("/home/users/qd201969/gtopo30_0.9x1.25.nc")
phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").load()
phis = phis/9.8 # transfer from m2/s2 to m

if cnlvl[0] < 0 :
    colr = cmaps.BlueDarkRed18(range(0,18,1))
    fcolors = colors.ListedColormap(np.vstack((colr[0:8,:],colr[10::,:])))
    print(fcolors.N)
else:
    fcolors = cmaps.precip2_17lev
cnlevels = np.arange(cnlvl[0], cnlvl[0]+cnlvl[1]*(fcolors.N-1), cnlvl[1])
norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=fcolors.N,extend='both')

fig = plt.figure(figsize=(12,12),dpi=100)
ax = fig.subplots(nrow, ncol, subplot_kw=dict(projection=ccrs.PlateCarree())) #sharex=True, sharey=True
for nl in range(0,len(lev),1):
    ds = xr.open_dataset('%s/ERA5_mon_vr_1979-2020.nc'%path)
    da = ds['vo'].sel(level=lev[nl],longitude=ilon,latitude=ilat,method="nearest").load()
    var = da.groupby(da.time.dt.month).mean('time').data*100000
    ds = xr.open_dataset('%s/ERA5_mon_z_1979-2020.nc'%path)
    da = ds['z'].sel(level=lev[nl],longitude=ilon,latitude=ilat,method="nearest").load()
    var1 = da.groupby(da.time.dt.month).mean('time').data/9.8
    '''
    ds = xr.open_dataset('%s/ERA5_mon_u_1979-2020.nc'%path)
    da = ds['u'].sel(level=lev[nl],longitude=ilon,latitude=ilat,method="nearest").load()
    uwnd = da.groupby(da.time.dt.month).mean('time').data
    ds = xr.open_dataset('%s/ERA5_mon_v_1979-2020.nc'%path)
    da = ds['v'].sel(level=lev[nl],longitude=ilon,latitude=ilat,method="nearest").load()
    vwnd = da.groupby(da.time.dt.month).mean('time').data
    var = dynamic_calc.calc_uv2vr_cfd(uwnd,vwnd,ilat,ilon) 
    '''
    del ds, da
    gc.collect()
    for nm in range(0,nrow,1):
        if nm == 0:
            shad = (var[0,:,:]+var[1,:,:]+var[11,:,:])/3.0
            cont1 = (var1[0,:,:]+var1[1,:,:]+var1[11,:,:])/3.0
        else:
            shad = np.mean(var[(3*nm-1):(3*nm+2),:,:],axis=0)
            cont1 = np.mean(var1[(3*nm-1):(3*nm+2),:,:],axis=0)
        axe = ax[nm][nl]
        axe.add_feature(cfeat.COASTLINE.with_scale('110m'),edgecolor='black', linewidth=0.8, zorder=1) 
        axe.set_title("%dhPa %s"%(lev[nl],titls[nm]),
            fontsize=title_font,fontdict=font)

        print('min:%f ; max:%f'%(np.nanmin(shad),np.nanmax(shad)))
        shad = axe.contourf(ilon, ilat, shad, cnlevels,
                     transform=ccrs.PlateCarree(),cmap=fcolors,extend='both',norm=norm)
        
        #wind = axe.quiver(ilon[::q_mis], ilat[::q_mis], uwnd1[::q_mis,::q_mis],vwnd1[::q_mis,::q_mis],
        #        pivot='mid',units='inches',scale=vcref[nl]*3,scale_units='inches',color="dimgray",
        #        width=0.02,headwidth=3,headlength=4.5,transform=ccrs.PlateCarree())

        cont = axe.contour(ilon, ilat, cont1, np.arange(1000,15000,cnlvl2[nl]), 
                     transform=ccrs.PlateCarree(), colors='darkviolet', linewidths=2)

        topo = axe.contour(ilon, ilat, phis, [1500,3000],
                     transform=ccrs.PlateCarree(),colors='black',linewidths=1.2)

        if nl == 0:
            axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
            axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
        if nm == (nrow-1):
            axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
            axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))

position = fig.add_axes([0.18, bmlo, 0.7, 0.01]) #left, bottom, width, height
cb = plt.colorbar(shad, cax=position ,orientation='horizontal')#, shrink=.9)
#axe.quiverkey(wind, 0.92, bmlo-0.01, vcref[nl], r'$%d m/s$'%vcref[nl], labelpos='N',coordinates='figure')

plt.figtext(0.02,bmlo-0.005, "%dhPa %s (%s)"%(lev[nl],drawvar[0],unit[0]), fontsize=title_font,
        horizontalalignment='left',verticalalignment='bottom')
plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
plt.savefig(figdir+"%s.png"%(drawvar[0]), bbox_inches='tight',pad_inches=0.01)

