'''
input start time, draw time and varname 
to draw two plots in one figure
left Njord d03, right PATH d04

20210925
'''

import xarray as xr
import numpy as np
import subprocess
import os 
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
from cartopy.io.shapereader import Reader
import cmaps
from PIL import Image, ImageDraw, ImageSequence
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, interplevel,ll_to_xy,xy_to_ll)

# constants
BIGFONT=22
MIDFONT=18
SMFONT=14

# define date of start time and forcast time
stime  = ['2021','09','15'] # year, month, date, hour
ftime  =[['2021','09','16','00'],
         ['2021','09','16','03'],
         ['2021','09','16','06'],
         ['2021','09','16','09']]
ftime0 =[['2021','09','15','21'], # use to calc accumulated pr,large than ftime
         ['2021','09','16','00'],
         ['2021','09','16','03'],
         ['2021','09','16','06']]

varname = ['slp','U10','V10']
drawvar = 'UV10 (m/s)'
#varname = ['RAINNC','RAINC']
#drawvar = 'preci(mm/3h)'
#varname = ['z']
lev = 500 #hPa

case = ['Njord','PATH']
path = ['/home/lzhenn/array74/Njord_Calypso/archive/',
        '/home/dataop/data/nmodel/wrf_fc/']
domin_select = ['d03','d04'] #njord, path
lat_sp = 1.0 #5.0 0.5#
lon_sp = 1.0 #10.00.5#
coast_shp = Reader(os.getenv('SHP_LIB')+'/china_coast/china_coastline.dbf').geometries()

# contour level for precip
#cnlevels = [0.1, 0.5, 1, 2, 3, 5, 8, 12, 16, 20, 25, 30, 40, 50, 70, 100, 150] #houly preci
#cnlevels = [0.1, 0.5, 1, 3, 5, 10,15,20, 30, 40, 60, 80, 100, 120, 150, 200, 250] #24h accumulated preci
#cnlevels = np.arange(1000,1017,1) #slp, hPa
#cnlevels = np.arange(586,603,1) #500Z, gpm
cnlevels = np.arange(0,8.5,0.5) #UV10 speeds
norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=cmaps.precip2_17lev.N,extend='both')
#norm = colors.Normalize(vmin=1000, vmax=1015)
NGITUDE_FORMATTER

ncfile = Dataset("/home/dataop/data/nmodel/wrf_fc/2021/202109/2021091512/wrfout_d04_2021-09-16_12:00:00")
llats  = getvar(ncfile,"XLAT")[0,0]
llatn  = getvar(ncfile,"XLAT")[-1,-1]
llonl  = getvar(ncfile,"XLONG")[0,0]
llonr  = getvar(ncfile,"XLONG")[-1,-1]

for t in range(len(ftime)):
    figtle = 'S'+'-'.join(stime)+' F'+'-'.join(ftime[t])+' '+drawvar
    figdir = '/home/lzhenn/cooperate/fig/S'+''.join(stime)+'.F'+''.join(ftime[t])+'.pr.png'
    fig = plt.figure(figsize=(12,9),dpi=100)
    for i in range(len(case)):
        if case[i] == case[0]: 
            domin = domin_select[0]
            q_mis=8 # wind vector plotting every q_miss grid
            filedir = path[i]+''.join(stime)+\
                    '/wrfout_'+domin+'_'+'-'.join(ftime[t][0:3])+'_'+ftime[t][3]+':00:00'
            filedir0= path[i]+''.join(stime)+\
                    '/wrfout_'+domin+'_'+'-'.join(ftime0[t][0:3])+'_'+ftime0[t][3]+':00:00'
        else:
            domin = domin_select[1]
            q_mis=8 # wind vector plotting every q_miss grid
            filedir = path[i]+stime[0]+'/'+''.join(stime[0:2])+'/'+''.join(stime)+\
                    '12/wrfout_'+domin+'_'+'-'.join(ftime[t][0:3])+'_'+ftime[t][3]+':00:00'
            filedir0= path[i]+stime[0]+'/'+''.join(stime[0:2])+'/'+''.join(stime)+\
                    '12/wrfout_'+domin+'_'+'-'.join(ftime0[t][0:3])+'_'+ftime0[t][3]+':00:00'

        # Open the NetCDF file
        ncfile = Dataset(filedir)
        x_y = to_np(ll_to_xy(ncfile,[llats,llatn],[llonl,llonr])) # return x and y
        var    = getvar(ncfile, varname[0])[x_y[1,0]:x_y[1,1],x_y[0,0]:x_y[0,1]] # first y, then x
        if len(var.dims) == 3:
            p      = getvar(ncfile, 'pressure')[x_y[1,0]:x_y[1,1],x_y[0,0]:x_y[0,1]] #hPa
            varnp  = to_np(interplevel(var, p, lev))/9.8
        else:
            varnp  = to_np(var)

        if varname[0] == 'RAINNC':
            print(''.join(ftime0[t])+' to '+''.join(ftime[t])+' preci')
            ncfile0 = Dataset(filedir0)
            varnp   = varnp+to_np(getvar(ncfile,'RAINC'))[x_y[1,0]:x_y[1,1],x_y[0,0]:x_y[0,1]]-\
                    to_np(getvar(ncfile0,"RAINC")[x_y[1,0]:x_y[1,1],x_y[0,0]:x_y[0,1]])-\
                    to_np(getvar(ncfile0,"RAINNC")[x_y[1,0]:x_y[1,1],x_y[0,0]:x_y[0,1]])

        # Get the latitude and longitude points
        lats, lons = latlon_coords(var)

        # Get the cartopy mapping object
        cart_proj = get_cartopy(var)

        # Set the GeoAxes to the projection used by WRF
        axe = plt.subplot(1,2,i+1,projection=cart_proj)    #创建子图
        
        coast_shp = Reader(os.getenv('SHP_LIB')+'/china_coast/china_coastline.dbf').geometries()
        coastline = cfeat.ShapelyFeature(coast_shp, ccrs.PlateCarree(), edgecolor='black', facecolor='none')
        axe.add_feature(coastline, linewidth=0.8,zorder=1)

        if len(varname) == 3:
            uvarnp = to_np(getvar(ncfile,varname[1]))[x_y[1,0]:x_y[1,1],x_y[0,0]:x_y[0,1]] # first y, then x
            vvarnp = to_np(getvar(ncfile,varname[2]))[x_y[1,0]:x_y[1,1],x_y[0,0]:x_y[0,1]] # first y, then x
            varnp = np.power((np.power(uvarnp,2)+np.power(vvarnp,2)),0.5) # speed of wind
            cont = axe.contourf(to_np(lons), to_np(lats), varnp, cnlevels, 
                         transform=ccrs.PlateCarree(),cmap=cmaps.precip2_17lev,extend='both',norm=norm)
            quv = axe.quiver(to_np(lons[::q_mis,::q_mis]),to_np(lats[::q_mis,::q_mis]),
                       uvarnp[::q_mis,::q_mis],vvarnp[::q_mis,::q_mis],zorder=2,
                       pivot='mid',units='inches',scale=30,scale_units='inches',color="dimgray",
                       width=0.02,headwidth=3,headlength=4.5,transform=ccrs.PlateCarree())

        else:
            cont = axe.contourf(to_np(lons), to_np(lats), varnp, cnlevels, 
                    transform=ccrs.PlateCarree(),cmap=cmaps.precip2_17lev,extend='both',norm=norm)

        if varname[0] == 'RAINNC':
            maxp = np.max(varnp)
            my,mx = np.where(varnp==maxp)
            print(my)
            mlon = to_np(lons[my,mx])
            mlat = to_np(lats[my,mx])
            #axe.plot(mlon,mlat,'ko',ms=10)
            axe.plot(lons[my,mx],lats[my,mx],'ko',ms=5,transform=ccrs.PlateCarree())
            axe.text(0.95,0.05, "Max pricip rate: %.2fmm/3hr @ (%.2fE,%.2fN)"%(maxp,mlon,mlat),
                    horizontalalignment='right',verticalalignment='bottom',transform=axe.transAxes)

        # Set the map bounds
        axe.set_xlim(cartopy_xlim(var))
        axe.set_ylim(cartopy_ylim(var))
        if case[i] == case[0]:
            axe.set_xticks(np.arange(np.ceil(lons[0,0]),np.ceil(lons[0,-1]),lat_sp), crs=ccrs.PlateCarree())
            axe.set_yticks(np.arange(np.ceil(lats[0,0]),np.ceil(lats[-1,0]),lon_sp), crs=ccrs.PlateCarree())
            axe.xaxis.set_major_formatter(LONGITUDE_FORMATTER) 
            axe.yaxis.set_major_formatter(LATITUDE_FORMATTER)
            axe.xaxis.set_major_formatter(LongitudeFormatter())
            axe.yaxis.set_major_formatter(LatitudeFormatter())
        
        axe.set_title(case[i],fontsize=MIDFONT) 

    fig.subplots_adjust(bottom=0.5,wspace=0.1,hspace=0.01)
    position = fig.add_axes([0.2, 0.45, 0.6, 0.015]) #left, bottom, width, height
    cb = plt.colorbar(cont, cax=position ,orientation='horizontal')#, shrink=.9)

    if len(varname) == 3:
        plt.quiverkey(quv, 0.87, 0.45, 10, r'$10 m/s$', labelpos='N',
                       coordinates='figure')

    plt.suptitle(figtle,x=0.5,y=0.95,fontsize=MIDFONT)
    plt.savefig(figdir,bbox_inches='tight',pad_inches=0.05)
    plt.show()

