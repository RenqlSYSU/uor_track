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
                 cartopy_ylim, latlon_coords, interplevel)

# constants
BIGFONT=22
MIDFONT=18
SMFONT=14

# define date of start time and forcast time
stime  = ['2021','06','22'] # year, month, date, hour
ftime  = ['2021','06','22','16']
ftime0 = ['2021','06','23','16'] # use to calc accumulated pr,large than ftime

#varname = ['slp','U10','V10']
#figtle = ''.join(ftime)+'_SLP(hPa)&UV10'
varname = ['RAINNC','RAINC']
drawvar = 'preci(mm/24h)'
#varname = ['z']
#figtle = ''.join(ftime)+'_500Z(gpm)'
lev = 500 #hPa

case = ['Njord','PATH']
domin_select = ['d03','d04'] #njord, path
lat_sp = 1.0 #5.0
lon_sp = 1.0 #10.0

path = '/home/lzhenn/cooperate/data/'
figtle = 'S'+''.join(stime)+' F'+''.join(ftime)+' '+drawvar
figdir = '/home/lzhenn/cooperate/fig/S'+''.join(stime)+'.F'+''.join(ftime)+'.pr.png'

coast_shp = Reader(os.getenv('SHP_LIB')+'/china_coast/china_coastline.dbf').geometries()

# contour level for precip
#cnlevels = [0.1, 0.5, 1, 2, 3, 5, 8, 12, 16, 20, 25, 30, 40, 50, 70, 100, 150] #houly preci
cnlevels = [0.1, 0.5, 1, 3, 5, 10,15,20, 30, 40, 60, 80, 100, 120, 150, 200, 250] #24h accumulated preci
#cnlevels = np.arange(1000,1017,1) #slp, hPa
#cnlevels = np.arange(586,603,1) #500Z, gpm
norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=cmaps.precip2_17lev.N,extend='both')
#norm = colors.Normalize(vmin=1000, vmax=1015)

fig = plt.figure(figsize=(12,9),dpi=100)
        
for i in range(len(case)):
    if case[i] == case[0]: 
        domin = domin_select[0]
        q_mis=15 # wind vector plotting every q_miss grid
        filedir = path+case[i]+'/'+''.join(stime)+\
                '/wrfout_'+domin+'_'+'-'.join(ftime[0:3])+'_'+ftime[3]+':00:00'
        filedir0= path+case[i]+'/'+''.join(stime)+\
                '/wrfout_'+domin+'_'+'-'.join(ftime0[0:3])+'_'+ftime0[3]+':00:00'
    else:
        domin = domin_select[1]
        q_mis=10 # wind vector plotting every q_miss grid
        filedir = path+case[i]+'/'+''.join(stime[0:2])+'/'+''.join(stime)+\
                '12/wrfout_'+domin+'_'+'-'.join(ftime[0:3])+'_'+ftime[3]+':00:00'
        filedir0= path+case[i]+'/'+''.join(stime[0:2])+'/'+''.join(stime)+\
                '12/wrfout_'+domin+'_'+'-'.join(ftime0[0:3])+'_'+ftime0[3]+':00:00'

    # Open the NetCDF file
    ncfile = Dataset(filedir)
    var    = getvar(ncfile, varname[0])
    if len(var.dims) == 3:
        p      = getvar(ncfile, 'pressure') #hPa
        varnp  = to_np(interplevel(var, p, lev))/9.8
    else:
        varnp  = to_np(var)

    if varname[0] == 'RAINNC':
        print(''.join(ftime)+' to '+''.join(ftime0)+' preci')
        ncfile0 = Dataset(filedir0)
        varnp   = to_np(getvar(ncfile0,"RAINC")) + to_np(getvar(ncfile0,"RAINNC"))\
                  -varnp-to_np(getvar(ncfile,'RAINC'))

    # Get the latitude and longitude points
    lats, lons = latlon_coords(var)

    # Get the cartopy mapping object
    cart_proj = get_cartopy(var)

    # Set the GeoAxes to the projection used by WRF
    axe = plt.subplot(1,2,i+1,projection=cart_proj)    #创建子图
    
    coast_shp = Reader(os.getenv('SHP_LIB')+'/china_coast/china_coastline.dbf').geometries()
    coastline = cfeat.ShapelyFeature(coast_shp, ccrs.PlateCarree(), edgecolor='grey', facecolor='none')
    axe.add_feature(coastline, linewidth=0.8)

    #cont = axe.contourf(to_np(lons), to_np(lats), varnp, 17,
    #             transform=ccrs.PlateCarree(),cmap=cmaps.precip2_17lev,extend='both')
    cont = axe.contourf(to_np(lons), to_np(lats), varnp, cnlevels, 
                 transform=ccrs.PlateCarree(),cmap=cmaps.precip2_17lev,extend='both',norm=norm)
    
    if len(varname) == 3:
        uvarnp = to_np(getvar(ncfile,varname[1]))
        vvarnp = to_np(getvar(ncfile,varname[2]))
        quv = axe.quiver(to_np(lons[::q_mis,::q_mis]),to_np(lats[::q_mis,::q_mis]),
                   uvarnp[::q_mis,::q_mis],vvarnp[::q_mis,::q_mis],
                   pivot='mid',units='inches',scale=30,scale_units='inches',
                   width=0.015,headwidth=3,headlength=4.5,transform=ccrs.PlateCarree())

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
    plt.quiverkey(quv, 0.87, 0.5, 10, r'$10 m/s$', labelpos='N',
                   coordinates='figure')

plt.suptitle(figtle,x=0.5,y=0.95,fontsize=MIDFONT)
plt.savefig(figdir,bbox_inches='tight',pad_inches=0.0)
plt.show()

