#!/usr/bin/env python
'''
input data path, varname to draw one figure

20210925
'''
import xarray as xr
import numpy as np
import subprocess
import os 
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
from cartopy.io.shapereader import Reader
import cmaps

# constants
BIGFONT=22
MIDFONT=18
SMFONT=14

varname = ["dtm"]
lat_sp = 0.1 #1.0 #
lon_sp = 0.1  #1.0 #
figdir = "/home/lzhenn/cooperate/fig/"+varname[0]+".png"
data  = ["/home/lzhenn/cooperate/data/Whole_HK_DTM_5m.nc",\
         "/home/lzhenn/cooperate/data/Whole_HK_DTM_100m.nc"]

# contour level for precip
cnlevels = np.arange(0,340,20) 
norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=cmaps.precip2_17lev.N,extend="both")

fig = plt.figure(figsize=(12,9),dpi=100)
# Open the NetCDF file
for nc in range(0,len(data),1):
    f  = xr.open_dataset(data[nc]) #yc,xc
    lat = f['lat'].data
    lon = f['lon'].data
    var  = f[varname[0]].data
    print(var)
    print(lat)

    term = data[nc].split("/")
    title=term[-1].split(".")

    axe = plt.subplot(1,2,nc+1,projection=ccrs.PlateCarree())    #创建子图
    axe.set_title(title[0],fontsize=MIDFONT) 

    coast_shp = Reader(os.getenv("SHP_LIB")+"/china_coast/china_coastline.dbf").geometries()
    coastline = cfeat.ShapelyFeature(coast_shp, ccrs.PlateCarree(), edgecolor="black", facecolor="none")
    axe.add_feature(coastline, linewidth=0.8,zorder=1)
    
    if nc == 1:
        cont = axe.pcolormesh(lon, lat, var, 
                transform=ccrs.PlateCarree(),cmap=cmaps.precip2_17lev,norm=norm)
    else:
        cont = axe.contourf(lon, lat, var, cnlevels, 
                     transform=ccrs.PlateCarree(),cmap=cmaps.precip2_17lev,extend='both',norm=norm)

    # Set the map bounds
    axe.set_xticks(np.arange(round(lon[0],1),round(lon[-1],1),lon_sp), crs=ccrs.PlateCarree())
    axe.set_yticks(np.arange(round(lat[-1],1),round(lat[0],1),lat_sp), crs=ccrs.PlateCarree())
    axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=""))
    axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=""))
    #axe.gridlines(draw_labels=True)
    #axe.xaxis.set_major_formatter(LONGITUDE_FORMATTER) 
    #axe.yaxis.set_major_formatter(LATITUDE_FORMATTER)

fig.subplots_adjust(bottom=0.5,wspace=0.15,hspace=0.01) # the wspace is 
#the width of the padding between subplots, as a fraction of the average Axes width
position = fig.add_axes([0.2, 0.45, 0.6, 0.015]) #left, bottom, width, height
cb = plt.colorbar(cont, cax=position ,orientation='horizontal')#, shrink=.9)
plt.savefig(figdir,bbox_inches="tight",pad_inches=0.1)


