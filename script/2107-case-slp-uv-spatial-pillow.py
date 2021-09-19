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
#stime  = ['2021','06','22','12'] # year, month, date, hour

#varname = ['slp','U10','V10']
#figtle = ''.join(ftime)+'_SLP(hPa)&UV10'
varname = ['RAINNC','RAINC']
drawvar = 'preci(mm/h)'
#varname = ['z']
#figtle = ''.join(ftime)+'_500Z(gpm)'
lev = 500 #hPa

domin_select = ['d03','d04'] #njord, path
lat_sp = 1.0 #5.0
lon_sp = 1.0 #10.0

coast_shp = Reader(os.getenv('SHP_LIB')+'/china_coast/china_coastline.dbf').geometries()

# contour level for precip
cnlevels = [0.1, 0.5, 1, 2, 3, 5, 8, 12, 16, 20, 25, 30, 40, 50, 70, 100, 150] #houly preci
#cnlevels = [0.1, 0.5, 1, 3, 5, 10,15,20, 30, 40, 60, 80, 100, 120, 150, 200, 250] #24h accumulated preci
#cnlevels = np.arange(1000,1017,1) #slp, hPa
#cnlevels = np.arange(1003,1011.5,0.5) #slp, hPa
#cnlevels = np.arange(586,603,1) #500Z, gpm
norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=cmaps.precip2_17lev.N,extend='both')
#norm = colors.Normalize(vmin=1000, vmax=1015)

path = '/home/lzhenn/cooperate/data/case_study/210628_blackrain/'
fn_stream = subprocess.check_output('ls '+path, shell=True).decode('utf-8')
case = fn_stream.split()
print(case)
print('case number : '+str(len(case)))

for nc in range(1,len(case),1):
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

    for itm in fn_list:
        fn_splt = itm.split('_')
        print(fn_splt) #['wrfout', 'd03', '2021-06-22', '11:00:00']
        figtle = case[nc]+'\n'+fn_splt[-2]+' '+fn_splt[-1]+' '+drawvar
        figdir = '/home/lzhenn/cooperate/fig/'+case[nc]+\
                 ''.join(fn_splt[-2:])+'.'+varname[0]+'.png'

        # Open the NetCDF file
        ncfile = Dataset(itm)
        var    = getvar(ncfile, varname[0])
        if len(var.dims) == 3:
            p      = getvar(ncfile, 'pressure') #hPa
            varnp  = to_np(interplevel(var, p, lev))/9.8
        else:
            varnp  = to_np(var)

        if varname[0] == 'RAINNC':
            if itm == fn_list[0]:
                varnp1 = to_np(getvar(ncfile,"RAINC")) + to_np(getvar(ncfile,"RAINNC"))
                varnp  = to_np(getvar(ncfile,"RAINC")) + to_np(getvar(ncfile,"RAINNC"))-varnp1
            else:
                varnp  = to_np(getvar(ncfile,"RAINC")) + to_np(getvar(ncfile,"RAINNC"))-varnp1
                varnp1 = to_np(getvar(ncfile,"RAINC")) + to_np(getvar(ncfile,"RAINNC"))

        # Get the latitude and longitude points
        lats, lons = latlon_coords(var)

        # Get the cartopy mapping object
        cart_proj = get_cartopy(var)

        # Set the GeoAxes to the projection used by WRF
        fig = plt.figure(figsize=(9,9),dpi=100)
        
        axe = plt.axes(projection=cart_proj)
        
        coast_shp = Reader(os.getenv('SHP_LIB')+'/china_coast/china_coastline.dbf').geometries()
        coastline = cfeat.ShapelyFeature(coast_shp, ccrs.PlateCarree(), edgecolor='grey', facecolor='none')
        axe.add_feature(coastline, linewidth=0.8)

        #cont = axe.contourf(to_np(lons), to_np(lats), varnp, 17,
        #             transform=ccrs.PlateCarree(),cmap=cmaps.precip2_17lev,extend='both')
        cont = axe.contourf(to_np(lons), to_np(lats), varnp, cnlevels, 
                     transform=ccrs.PlateCarree(),cmap=cmaps.precip2_17lev,extend='both',norm=norm)
        plt.colorbar(cont, ax=axe, shrink=.7)
        
        if len(varname) == 3:
            uvarnp = to_np(getvar(ncfile,varname[1]))
            vvarnp = to_np(getvar(ncfile,varname[2]))
            quv = axe.quiver(to_np(lons[::q_mis,::q_mis]),to_np(lats[::q_mis,::q_mis]),
                       uvarnp[::q_mis,::q_mis],vvarnp[::q_mis,::q_mis],
                       pivot='mid',units='inches',scale=30,scale_units='inches',
                       width=0.015,headwidth=3,headlength=4.5,transform=ccrs.PlateCarree())
            plt.quiverkey(quv, 1.08, -0.05, 10, r'$10 m/s$', labelpos='N',
                           coordinates='axes')

        # Set the map bounds
        axe.set_xlim(cartopy_xlim(var))
        axe.set_ylim(cartopy_ylim(var))
        if itmsplt[0] == 'njord':
            axe.set_xticks(np.arange(np.ceil(lons[0,0]),np.floor(lons[0,-1]),lat_sp), crs=ccrs.PlateCarree())
            axe.set_yticks(np.arange(np.ceil(lats[0,0]),np.floor(lats[-1,0]),lon_sp), crs=ccrs.PlateCarree())
            axe.xaxis.set_major_formatter(LONGITUDE_FORMATTER) 
            axe.yaxis.set_major_formatter(LATITUDE_FORMATTER)
            axe.xaxis.set_major_formatter(LongitudeFormatter())
            axe.yaxis.set_major_formatter(LatitudeFormatter())
        
        plt.title(figtle,fontsize=MIDFONT)
        plt.savefig(figdir,bbox_inches='tight',pad_inches=0.0)
        plt.show()
        del fig, axe

    fn_stream = subprocess.check_output('ls /home/lzhenn/cooperate/fig/'+\
                 case[nc]+'*.'+varname[0]+'.png', shell=True).decode('utf-8')
    fn_list   = fn_stream.split()
    print(fn_list[0])
    print('filenumber : '+str(len(fn_list)))
    gif_name = '/home/lzhenn/cooperate/fig/'+case[nc]+'.'+varname[0]+'.gif'

    frames = []
    for itm in fn_list:
        frame = Image.open(itm)
        frames.append(frame)

    frames[0].save(gif_name, save_all=True, append_images=frames[1:],\
                duration = 1000, loop=0, disposal=1)
    subprocess.run('rm -f /home/lzhenn/cooperate/fig/'+case[nc]+'*.'+varname[0]+'.png',shell=True)


