'''
input data path, varname to draw one figure

20210925
'''
import xarray as xr
import numpy as np
import subprocess
import os 
from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
from cartopy.io.shapereader import Reader
import cmaps
from datetime import datetime
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, interplevel)
matplotlib.use('Agg')

# constants
BIGFONT=22
MIDFONT=18
SMFONT=14

# define date of start time and forcast time
ftime  = datetime(2018,9,14,00)
ftime0 = datetime(2018,9,17,00)
domain = 'd02'
#data = "/home/metctm1/array_hq133/data/s2s/wrfonly/clim/2011070100/wrfout_d01_2011-07-23_00:00:00"
path = '/home/lzhenn/cooperate/data/case_study/coupled/2018091200'
data = "%s/wrfout_%s_%s:00:00"%(path,domain,ftime.strftime('%Y-%m-%d_%H')) 

#varname = ['slp','U10','V10']
#drawvar = 'UV10 (m/s)'
varname = ['RAINNC','RAINC']
drawvar = 'preci(mm)'
#varname = ['z']
#drawvar = '500 Z'
figtle = 'CTRL'
lev = 500 #hPa
figdir = '/home/lzhenn/cooperate/fig/%s_%s.png'%(ftime.strftime('%Y-%m-%d_%H'),varname[0])

lat_sp = 1.0 #5.0
lon_sp = 1.0 #10.0
q_mis=15 # wind vector plotting every q_miss grid

# contour level for precip
#cnlevels = [0.1, 0.5, 1, 2, 3, 5, 8, 12, 16, 20, 25, 30, 40, 50, 70, 100, 150] #houly preci
cnlevels = [0.1, 0.5, 1, 3, 5, 10,15,20, 30, 40, 60, 80, 100, 120, 150, 200, 250] #24h accumulated preci
#cnlevels = np.arange(1000,1017,1) #slp, hPa
#cnlevels = np.arange(586,603,1) #500Z, gpm
#cnlevels = np.arange(0,17,1) #UV10 speeds
norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=cmaps.precip2_17lev.N,extend='both')
#norm = colors.Normalize(vmin=1000, vmax=1015)

# Open the NetCDF file
ncfile = Dataset(data)
var    = getvar(ncfile, varname[0])
if len(var.dims) == 3:
    p      = getvar(ncfile, 'pressure') #hPa
    varnp  = to_np(interplevel(var, p, lev))/9.8
else:
    varnp  = to_np(var)

if varname[0] == 'RAINNC':
    filedir0 = "%s/wrfout_%s_%s:00:00"%(path,domain,ftime0.strftime('%Y-%m-%d_%H')) 
    print('%s to %s preci'%(ftime.strftime('%Y-%m-%d_%H'),ftime0.strftime('%Y-%m-%d_%H')))
    ncfile0 = Dataset(filedir0)
    varnp   = to_np(getvar(ncfile0,"RAINC")) + to_np(getvar(ncfile0,"RAINNC"))\
              -varnp-to_np(getvar(ncfile,'RAINC'))

# Get the latitude and longitude points
lats, lons = latlon_coords(var)

# Get the cartopy mapping object
cart_proj = get_cartopy(var)
print(cart_proj)

# Set the GeoAxes to the projection used by WRF
fig = plt.figure(figsize=(12,9),dpi=100)
axe = plt.axes(projection=cart_proj)
#axe = plt.subplot(1,2,i+1,projection=cart_proj)    #创建子图

#coast_shp = Reader(os.getenv('SHP_LIB')+'/china_coast/china_coastline.dbf').geometries()
#coastline = cfeat.ShapelyFeature(coast_shp, ccrs.PlateCarree(), edgecolor='black', facecolor='none')
#axe.add_feature(coastline, linewidth=0.8,zorder=1)
province_shp_file=os.getenv('SHP_LIB')+'/cnmap/cnhimap.dbf'
province_shp=shpreader.Reader(province_shp_file).geometries()
axe.add_geometries(province_shp, ccrs.PlateCarree(),facecolor='none', edgecolor='black',linewidth=1., zorder = 1)
if len(varname) == 3:
    uvarnp = to_np(getvar(ncfile,varname[1]))
    vvarnp = to_np(getvar(ncfile,varname[2]))
    varnp = np.power((np.power(uvarnp,2)+np.power(vvarnp,2)),0.5) # speed of wind
    cont = axe.pcolormesh(to_np(lons), to_np(lats), varnp, 
            transform=ccrs.PlateCarree(),cmap=cmaps.precip2_17lev,norm=norm)
    quv = axe.quiver(to_np(lons[::q_mis,::q_mis]),to_np(lats[::q_mis,::q_mis]),
               uvarnp[::q_mis,::q_mis],vvarnp[::q_mis,::q_mis],zorder=2,
               pivot='mid',units='inches',scale=30,scale_units='inches',color="dimgray",
               width=0.02,headwidth=3,headlength=4.5,transform=ccrs.PlateCarree())
    axe.quiverkey(quv, 1.05, 0.97, 10, r'$10 m/s$', labelpos='N',
               coordinates='axes')
else:
    cont = axe.contourf(to_np(lons), to_np(lats), varnp, cnlevels, 
            transform=ccrs.PlateCarree(),cmap=cmaps.precip2_17lev,extend='both',norm=norm)

# Set the map bounds
axe.set_xlim(cartopy_xlim(var))
axe.set_ylim(cartopy_ylim(var))
#axe.gridlines(draw_labels=True)
#axe.set_xticks(np.arange(np.ceil(lons[0,0]),np.ceil(lons[0,-1]),lat_sp), crs=ccrs.PlateCarree())
#axe.set_yticks(np.arange(np.ceil(lats[0,0]),np.ceil(lats[-1,0]),lon_sp), crs=ccrs.PlateCarree())
#axe.xaxis.set_major_formatter(LONGITUDE_FORMATTER) 
#axe.yaxis.set_major_formatter(LATITUDE_FORMATTER)
#axe.xaxis.set_major_formatter(LongitudeFormatter())
#axe.yaxis.set_major_formatter(LatitudeFormatter())

axe.set_title(figtle,fontsize=MIDFONT) 
plt.colorbar(cont, ax=axe, shrink=.7)
plt.savefig(figdir,bbox_inches='tight',pad_inches=0.1)
plt.show()

