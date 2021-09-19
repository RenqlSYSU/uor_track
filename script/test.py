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
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                         cartopy_ylim, latlon_coords)

# constants
BIGFONT=22
MIDFONT=18
SMFONT=14

# define date of start time and forcast time
stime  = ['2021','06','22'] # year, month, date, hour
ftime1 = ['2021','06','23','08']
ftime2 = ['2021','06','24','08']

path = '/home/lzhenn/cooperate/data'
filedir1 = ['/Njord/'+''.join(stime)+'/wrfout_d03_'+'-'.join(ftime1[0:3])+'_'+ftime1[3]+':00:00',\
            '/PATH/'+''.join(stime[0:2])+'/'+''.join(stime)+'12/wrfout_d04_'+'-'.join(ftime1[0:3])+'_'+ftime1[3]+':00:00']
filedir2 = ['/Njord/'+''.join(stime)+'/wrfout_d03_'+'-'.join(ftime2[0:3])+'_'+ftime2[3]+':00:00',\
            '/PATH/'+''.join(stime[0:2])+'/'+''.join(stime)+'12/wrfout_d04_'+'-'.join(ftime2[0:3])+'_'+ftime2[3]+':00:00']
figtle = ['Njord','PATH','S'+''.join(stime)+'.F'+''.join(ftime1)+' 24h Accumulated Precipitation (mm)']
figdir = '/home/lzhenn/cooperate/script/S'+''.join(stime)+'.F'+''.join(ftime1)+'.pr.png'

coast_shp = Reader(os.getenv('SHP_LIB')+'/china_coast/china_coastline.shp').geometries()

cnlevels = [0.1, 0.5, 1, 2, 3, 5, 8, 12, 16, 20, 25, 30, 40, 50, 70, 100, 150]
norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=cmaps.precip2_17lev.N,extend='both')

# Create a figure
fig = plt.figure(figsize=(12,6))

for i in range(2):
    # Open the NetCDF file
    ncfile = Dataset(path+filedir1[i])
    var    = getvar(ncfile, "RAINC" )
    nproj  = str(var.attrs['projection']).split('(')
    print(nproj)

    # Get the latitude and longitude points
    lats, lons = latlon_coords(var)

    # Get the cartopy mapping object
    cart_proj = get_cartopy(var)

    # Set the GeoAxes to the projection used by WRF
    #axe = plt.axes(projection=cart_proj)

    ncfile = Dataset(path+filedir2[i])
    varnp  = to_np(getvar(ncfile,"RAINC")) + to_np(getvar(ncfile,"RAINNC"))
    ncfile = Dataset(path+filedir1[i])
    varnp  = varnp - (to_np(getvar(ncfile,"RAINC")) + to_np(getvar(ncfile,"RAINNC")))

    axe = plt.subplot(1,2,i+1,projection=cart_proj)    #创建子图
    # Download and add the states and coastlines
    #states = cfeat.NaturalEarthFeature(category="cultural", scale="10m",facecolor="none", name="admin_1_states_provinces_shp")
    #axe.add_feature(states, linewidth=.5, edgecolor="black")
    #axe.coastlines('50m', linewidth=0.8)
    #axe.add_feature(cfeat.COASTLINE.with_scale('10m'), linewidth=1,color='k')
    axe.add_geometries(coast_shp, ccrs.PlateCarree(),facecolor='none', edgecolor='grey',linewidth=.5, zorder = 1)
    #coastline = cfeat.ShapelyFeature(coast_shp, ccrs.PlateCarree(), edgecolor='grey', facecolor='none')
    #axe.add_feature(coastline, linewidth=0.5)
    #LAKES_border = cfeat.NaturalEarthFeature('physical', 'lakes', '10m', edgecolor='k', facecolor='never')
    #axe.add_feature(LAKES_border, linewidth=0.8)
    #axe.add_feature(cfeat.LAND.with_scale('10m'))
    #axe.add_feature(cfeat.OCEAN.with_scale('10m'))
    #axe.add_feature(cfeat.RIVERS.with_scale('10m'))

    # Make the contour outlines and filled contours for variable
    #plt.contour(to_np(lons), to_np(lats), to_np(var), 10, colors="black",transform=ccrs.PlateCarree())
    #plt.contourf(to_np(lons), to_np(lats), varnp, 17, transform=ccrs.PlateCarree(),cmap=cmaps.precip2_17lev)
    plt.contourf(to_np(lons), to_np(lats), varnp, cnlevels, 
                 transform=ccrs.PlateCarree(),cmap=cmaps.precip2_17lev,norm=norm,extend='both')

    # Add a color bar
    plt.colorbar(ax=axe, shrink=.58)
    #plt.colorbar(ax=axe, drawedges=True, orientation='vertical',spacing='uniform')

    # Set the map bounds
    axe.set_xlim(cartopy_xlim(var))
    axe.set_ylim(cartopy_ylim(var))
    #axe.set_xticks(np.arange(np.ceil(lons[0,0]),np.floor(lons[0,-1]),1.0), crs=ccrs.PlateCarree())
    #axe.set_yticks(np.arange(np.ceil(lats[0,0]),np.floor(lats[-1,0]),1.0), crs=ccrs.PlateCarree())
    
    # Define gridline locations and draw the lines using cartopy's built-in gridliner:
    xticks = list(np.arange(np.ceil(lons[0,0]),np.floor(lons[0,-1]),1.0))
    yticks = list(np.arange(np.ceil(lats[0,0]),np.floor(lats[-1,0]),1.0))
    axe.gridlines(xlocs=xticks, ylocs=yticks)

    # Label the end-points of the gridlines using the custom tick makers:
    axe.xaxis.set_major_formatter(LONGITUDE_FORMATTER) 
    axe.yaxis.set_major_formatter(LATITUDE_FORMATTER)
   
    if nproj[0] == 'LambertConformal':
        # xticks and yticks only is list
        Lambert_tick.lambert_xticks(axe, xticks)  
        Lambert_tick.lambert_yticks(axe, yticks)
    else:
        axe.xaxis.set_major_formatter(LongitudeFormatter())
        axe.yaxis.set_major_formatter(LatitudeFormatter())

    # Add the gridlines
    #axe.gridlines(color="black", linestyle="dotted")

    plt.title(figtle[i],fontsize=SMFONT)

plt.suptitle(figtle[2],fontsize=MIDFONT)
plt.savefig(figdir)
plt.show()

