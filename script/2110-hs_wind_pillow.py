#!/usr/bin/env python
'''
input data path, varname to draw one figure

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

# constants
BIGFONT=22
MIDFONT=18
SMFONT=14

data  = "/home/lzhenn/Calypso/GBA_hsig_d02.nc"
data1 = "/home/lzhenn/Calypso/GBA_wind_d02.nc"

varname = ["hs","xwnd","ywnd"]
drawvar = "sea surface wave significant height (cm) & surface wind"
#varname = ["slp","U10","V10"]
#drawvar = "UV10 (m/s)"
#varname = ["RAINNC","RAINC"]
#drawvar = "preci(mm/24h)"
#varname = ["z"]
#figtle = "".join(ftime)+"_500Z(gpm)"
lev = 500 #hPa
lonl= 112.5 #0  #
lonr= 115   #360#
lats= 21.5  #0  #
latn= 23    #90 #
lat_sp = 0.5 #1.0 #
lon_sp = 0.5  #1.0 #
q_mis=15 # wind vector plotting every q_miss grid

# contour level for precip
#cnlevels = [0.1, 0.5, 1, 2, 3, 5, 8, 12, 16, 20, 25, 30, 40, 50, 70, 100, 150] #houly preci
#cnlevels = [0.1, 0.5, 1, 3, 5, 10,15,20, 30, 40, 60, 80, 100, 120, 150, 200, 250] #24h accumulated preci
#cnlevels = np.arange(1000,1017,1) #slp, hPa
#cnlevels = np.arange(586,603,1) #500Z, gpm
cnlevels = np.arange(0,85,5) #UV10 speeds
norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=cmaps.precip2_17lev.N,extend="both")
#norm = colors.Normalize(vmin=1000, vmax=1015)

#=====================================================
# read data 
#====================================================================
# Open the NetCDF file
f  = xr.open_dataset(data) #yc,xc
lat = f['latitude'].sel(xc=0)
lon = f['longitude'].sel(yc=0)
ixc= f.xc[(lon>=lonl) & (lon<=lonr)]
iyc = f.yc[(lat>=lats) & (lat<=latn)]
time = f.time
var  = f[varname[0]].sel(yc=iyc,xc=ixc).load()
ilat = f['latitude'].sel(yc=iyc,xc=ixc).load()
ilon = f['longitude'].sel(yc=iyc,xc=ixc).load()
var=var*100
print(var)

f  = xr.open_dataset(data1) 
uvar  = f[varname[1]].sel(yc=iyc,xc=ixc).load()
vvar  = f[varname[2]].sel(yc=iyc,xc=ixc).load()

#=====================================================
# draw figure 
#====================================================================
for nt in range(0,len(time),1):
    figdir = "/home/lzhenn/cooperate/fig/"+varname[0]+str(nt)+".png"
    print(figdir)
    
    fig = plt.figure(figsize=(12,9),dpi=100)
    axe = plt.axes(projection=ccrs.PlateCarree())
    axe.set_title(drawvar+"\n"+str(time[nt].data).strip("0").strip("."),fontsize=MIDFONT) 

    coast_shp = Reader(os.getenv("SHP_LIB")+"/china_coast/china_coastline.dbf").geometries()
    coastline = cfeat.ShapelyFeature(coast_shp, ccrs.PlateCarree(), edgecolor="black", facecolor="none")
    axe.add_feature(coastline, linewidth=0.8,zorder=1)

    cont = axe.contourf(ilon, ilat, var[nt,:,:], cnlevels, 
            transform=ccrs.PlateCarree(),cmap=cmaps.precip2_17lev,extend="both",norm=norm)
    #cont = axe.pcolormesh(ilon, ilat, var, 
    #        transform=ccrs.PlateCarree(),cmap=cmaps.precip2_17lev,norm=norm)
    quv = axe.quiver(ilon[::q_mis,::q_mis],ilat[::q_mis,::q_mis],
               uvar[nt,::q_mis,::q_mis],vvar[nt,::q_mis,::q_mis],zorder=2,
               pivot="mid",units="inches",scale=30,scale_units="inches",color="dimgray",
               width=0.02,headwidth=3,headlength=4.5,transform=ccrs.PlateCarree())
    axe.quiverkey(quv, 1.05, 0.97, 10, r"$10 m/s$", labelpos="N",
               coordinates="axes")

    # Set the map bounds
    axe.set_xticks(np.arange(lonl,lonr+lon_sp,lon_sp), crs=ccrs.PlateCarree())
    axe.set_yticks(np.arange(lats,latn+lat_sp,lat_sp), crs=ccrs.PlateCarree())
    axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=""))
    axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=""))
    #axe.gridlines(draw_labels=True)
    #axe.xaxis.set_major_formatter(LONGITUDE_FORMATTER) 
    #axe.yaxis.set_major_formatter(LATITUDE_FORMATTER)

    plt.colorbar(cont, ax=axe, shrink=.55)
    plt.savefig(figdir,bbox_inches="tight",pad_inches=0.1)
    del fig, axe

#=====================================================
# animation
#====================================================================
fn_stream = subprocess.check_output("ls -rt /home/lzhenn/cooperate/fig/"+\
             varname[0]+"*.png", shell=True).decode("utf-8")
fn_list   = fn_stream.split()
print(fn_list[0])
print("filenumber : "+str(len(fn_list)))
gif_name = "/home/lzhenn/cooperate/fig/"+varname[0]+".gif"

frames = []
for itm in fn_list:
    frame = Image.open(itm)
    frames.append(frame)

frames[0].save(gif_name, save_all=True, append_images=frames[1:],\
            duration = 250, loop=0, disposal=0)
# duration: The time to display the current frame of the GIF, in milliseconds.
# loop: The number of times the GIF should loop. 0 means that it will loop forever.
subprocess.run("rm -f /home/lzhenn/cooperate/fig/"+varname[0]+"*.png",shell=True)


