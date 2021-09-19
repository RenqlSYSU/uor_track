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
                 cartopy_ylim, latlon_coords, interplevel)

# constants
BIGFONT=22
MIDFONT=18
SMFONT=14

# define date of start time and forcast time
stime  = ['2021','06','27','12'] # year, month, date, hour
ftime  = ['2021','06','27','12']
ftime0 = ['2021','06','28','08'] # use to calc accumulated pr,large than ftime

varname = ['slp'] #,'U10','V10'
drawvar = 'SLP(hPa)' #UV10
#varname = ['RAINNC','RAINC']
#drawvar = 'preci(mm/24h)'
#varname = ['z']
#drawvar = '500Z(gpm)'
lev = 500 #hPa
#figtle = ''.join(ftime0)+'_'+''.join(ftime)+'_'+drawvar

domin_select = ['d03','d04'] #njord, path
lat_sp = 1.0 #10.0 #
lon_sp = 1.0 # 5.0 #

path = '/home/lzhenn/cooperate/data/case_study/210628_blackrain/'
#fn_stream = subprocess.check_output('ls '+path, shell=True).decode('utf-8')
#case = fn_stream.split()
#print(case)
#print('case number : '+str(len(case)))
case = ['njord.gfs.062712.0d25','njord.gfs.062712.1d00','path.gfs.062712.1d00']
#case = ['njord.gfs.062712.0d25']

coast_shp = Reader(os.getenv('SHP_LIB')+'/china_coast/china_coastline.dbf').geometries()

# contour level for precip
#cnlevels = [0.1, 0.5, 1, 2, 3, 5, 8, 12, 16, 20, 25, 30, 40, 50, 70, 100, 150] #houly preci
#cnlevels = [0.1, 0.5, 1, 3, 5, 10,15,20, 30, 40, 60, 80, 100, 120, 150, 200, 250] #24h accumulated preci
#cnlevels = np.arange(1000,1017,1) #slp, hPa, d01
cnlevels = np.arange(1004,1007.4,0.2) #slp, hPa, d03
#cnlevels = np.arange(586,603,1) #500Z, gpm
norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=cmaps.precip2_17lev.N,extend='both')
#norm = colors.Normalize(vmin=1000, vmax=1015)

for nc in range(len(case)):
    itmsplt = case[nc].split('.')
    if itmsplt[0] == 'njord':
        domin = domin_select[0]
        q_mis=15 # wind vector plotting every q_miss grid
    else:
        domin = domin_select[1]
        q_mis=10 # wind vector plotting every q_miss grid

    fn_stream = subprocess.check_output('ls '+path+case[nc]+'/wrfout_'+domin+'_*', shell=True).decode('utf-8')
    fn_list   = fn_stream.split()
    print(fn_list[0])
    print('filenumber : '+str(len(fn_list)))
    
    figtle = case[nc]+' '+drawvar
    figdir = '/home/lzhenn/cooperate/fig/case_study-'+case[nc]+'.'+varname[0]+'.jpg'
    fig = plt.figure(figsize=(12,9),dpi=300)
        
    i = -1
    for itm in fn_list[0:6]:
        i = i+1
        fn_splt = itm.split('_')
        print(fn_splt) #['wrfout', 'd03', '2021-06-22', '11:00:00']

        # Open the NetCDF file
        ncfile = Dataset(itm)
        var    = getvar(ncfile, varname[0])
        if len(var.dims) == 3:
            p      = getvar(ncfile, 'pressure') #hPa
            varnp  = to_np(interplevel(var, p, lev))/9.8
        else:
            varnp  = to_np(var)

        if varname[0] == 'RAINNC':
            print(''.join(ftime0)+' to '+''.join(ftime)+' preci')
            filedir0 = path+case[nc]+'/wrfout_'+domin+'_'+'-'.join(ftime0[0:3])+'_'+ftime0[3]+':00:00'
            ncfile0 = Dataset(filedir0)
            varnp   = varnp + to_np(getvar(ncfile,'RAINC'))-\
                      (to_np(getvar(ncfile0,"RAINC")) + to_np(getvar(ncfile0,"RAINNC")))

        # Get the latitude and longitude points
        lats, lons = latlon_coords(var)

        # Get the cartopy mapping object
        cart_proj = get_cartopy(var)

        # Set the GeoAxes to the projection used by WRF
        axe = plt.subplot(2,3,i+1,projection=cart_proj)    #创建子图
        
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
        if itmsplt[0] == 'njord':
            axe.set_xticks(np.arange(np.ceil(lons[0,0]),np.ceil(lons[0,-1]),lat_sp), crs=ccrs.PlateCarree())
            axe.set_yticks(np.arange(np.ceil(lats[0,0]),np.ceil(lats[-1,0]),lon_sp), crs=ccrs.PlateCarree())
            axe.xaxis.set_major_formatter(LONGITUDE_FORMATTER) 
            axe.yaxis.set_major_formatter(LATITUDE_FORMATTER)
            axe.xaxis.set_major_formatter(LongitudeFormatter())
            axe.yaxis.set_major_formatter(LatitudeFormatter())
        
        axe.set_title(fn_splt[-2]+' '+fn_splt[-1],fontsize=SMFONT) 

    fig.subplots_adjust(bottom=0.1,wspace=0.2,hspace=0.01)
    position = fig.add_axes([0.15, 0.08, 0.7, 0.015]) #left, bottom, width, height
    cb = plt.colorbar(cont, cax=position ,orientation='horizontal')#, shrink=.9)
    if len(varname) == 3:
        plt.quiverkey(quv, 0.87, 0.08, 10, r'$10 m/s$', labelpos='N',
                       coordinates='figure')

    plt.suptitle(figtle,x=0.5,y=0.9,fontsize=MIDFONT)
    plt.savefig(figdir,bbox_inches='tight',pad_inches=0.0)
    plt.show()

