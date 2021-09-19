'''
try use ArtistAnimation to produce animation of hourly preci
but it does not work.
'''
import xarray as xr
import numpy as np
import subprocess
import os 
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
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

figtle = "Precipitation (mm/hour)"
case = ['Njord','PATH']
domin = ['03','04']
nc = 1

path = '/home/lzhenn/cooperate/data'
filedir= [path+'/Njord/'+''.join(stime)+'/',\
          path+'/PATH/'+''.join(stime[0:2])+'/'+''.join(stime)+'12/']
coast_shp = Reader(os.getenv('SHP_LIB')+'/china_coast/china_coastline.shp').geometries()

# contour level for precip
cnlevels = [0.1, 0.5, 1, 2, 3, 5, 8, 12, 16, 20, 25, 30, 40, 50, 70, 100, 150]
norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=cmaps.precip2_17lev.N,extend='both')

fn_stream = subprocess.check_output('ls '+filedir[nc]+'wrfout_d'+domin[nc]+'_*', shell=True).decode('utf-8')
fn_list   = fn_stream.split()
print(fn_list[0])
print('filenumber : '+str(len(fn_list)))

ncfile = Dataset(fn_list[0])
var    = getvar(ncfile, "RAINC")
varnp1 = to_np(getvar(ncfile,"RAINC")) + to_np(getvar(ncfile,"RAINNC"))
print(var)

# Get the latitude and longitude points
lats, lons = latlon_coords(var)

# Get the cartopy mapping object
cart_proj = get_cartopy(var)

fig = plt.figure(figsize=(12,6))
    
# Set the GeoAxes to the projection used by WRF
axe = plt.axes(projection=cart_proj)
axe.add_geometries(coast_shp, ccrs.PlateCarree(),facecolor='none', edgecolor='grey',linewidth=.5, zorder = 1)

artists = []
for itm in fn_list[1:10]:
    fn_splt = itm.split('_')
    print(fn_splt) #['wrfout', 'd03', '2021-06-22', '11:00:00']
    figtle = case[nc]+' '+fn_splt[2]+' '+fn_splt[3]+' precip(mm/h)'
    figdir = '/home/lzhenn/cooperate/fig/'+case[nc]+'.S'+''.join(stime)+ \
             '.F'+''.join(fn_splt[2].split('-'))+'-'+fn_splt[3]+'.pr.png'

    # Open the NetCDF file
    ncfile = Dataset(itm)
    varnp  = to_np(getvar(ncfile,"RAINC")) + to_np(getvar(ncfile,"RAINNC"))-varnp1
    varnp1 = to_np(getvar(ncfile,"RAINC")) + to_np(getvar(ncfile,"RAINNC"))

    plot = plt.contourf(to_np(lons), to_np(lats), varnp, cnlevels, 
                  transform=ccrs.PlateCarree(),cmap=cmaps.precip2_17lev,norm=norm,extend='both')

    plt.colorbar(ax=axe, shrink=.58)

    # Set the map bounds
    axe.set_xlim(cartopy_xlim(var))
    axe.set_ylim(cartopy_ylim(var))
    #axe.set_xticks(np.arange(np.ceil(lons[0,0]),np.floor(lons[0,-1]),1.0), crs=ccrs.PlateCarree())
    #axe.set_yticks(np.arange(np.ceil(lats[0,0]),np.floor(lats[-1,0]),1.0), crs=ccrs.PlateCarree())
    axe.xaxis.set_major_formatter(LongitudeFormatter())
    axe.yaxis.set_major_formatter(LatitudeFormatter())
    
    text = plt.title(figtle,fontsize=SMFONT)
    artists.append(plot)
    
ani = animation.ArtistAnimation(fig=fig, artists=artists, repeat=True, interval=100)
plt.show()
ani.save('./2.mp4', fps=30)
#plt.savefig(figdir)
#plt.show()

