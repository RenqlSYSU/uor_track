import xarray as xr
import pandas as pd
import numpy as np
import gc #garbage collector
import subprocess, os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import colors
from matplotlib import cm
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.io.shapereader import Reader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmaps, wrf
from netCDF4 import Dataset
from datetime import datetime
from PIL import Image

title_font=14
label_font=10
plt.rcParams["font.weight"] = "bold"
font = {'family': 'sans-serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black',
        }

figdir = "/home/lzhenn/cooperate/fig"
def main_run():
    '''
    path = "/home/metctm1/array_hq133/data/s2s/wrfroms/realtime/aegir_2022070100"
    filname = ['%s/wrfout_d01_2022-07-18_09:00:00'%path,
        '%s/njord_his_d01.20220725.nc'%path]
    draw_domain_terrain(filname,'WRF+ROMS Domain')
    '''
    
    paths = ['/home/metctm1/array_hq133/data/s2s/wrfonly/realtime',
             '/home/metctm1/array_hq133/data/s2s/wrfroms/realtime']
    indir = ['%s/WRF_free'%paths[0],
             '%s/CP_free'%paths[1]]
    case = [x.split('/')[-1] for x in indir] 
    ftime  = pd.date_range(start='2022-07-05 00',end='2022-07-10 23',
        freq='1H',closed=None)
    
    cnlevels = np.arange(20,31.6,0.1)
    ncmap = colors.ListedColormap(cmaps.rainbow(
        np.linspace(0,1,len(cnlevels)+1))) 
    for nc in range(1,len(case)):
        for nt in ftime:
            filname = '%s/wrfout_d01_%s'%(indir[nc],nt.strftime("%Y-%m-%d_%H:00:00"))
            draw_wrf_var(case[nc],filname,'SST','degC',cnlevels,ncmap)
        animat('%s/%s_%s_*.png'%(figdir,'SST',case[nc]),
            '%s/%s_%s.gif'%(figdir,'SST',case[nc]))

def draw_domain_terrain(filname,figtle):
    figname = '%s/domain.png'%(figdir)
    cnlevels = np.array([0, 20, 50, 100, 150, 200, 250, 300, 350, 
        400, 450, 500, 550, 600, 700, 800, 900, 1000, 1100, 1200, 
        1300, 1400, 1500, 1700, 2000, 2500, 3000, 3500, 4000, 4500])
    ncmap = cmaps.OceanLakeLandSnow(np.linspace(0, 1, len(cnlevels)+1))
    ncmap[0,:] = np.array([0.5,0.8,1.0,1.0])
    ncmap = colors.ListedColormap(ncmap)
    norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=ncmap.N, extend='both')
    
    d01 = Dataset(filname[0])
    hgt = wrf.getvar(d01,'HGT')
    lake = wrf.getvar(d01,'LAKEMASK').data #1 FOR LAKE, 0 FOR NON-LAKE
    land = wrf.getvar(d01,'LANDMASK').data #1 FOR LAND, 0 FOR WATER
    hgt.data = np.where(lake==1, -1, hgt.data)
    hgt.data = np.where(land==0, -2, hgt.data)
    del lake, land
    gc.collect()

    lats, lons = wrf.latlon_coords(hgt)
    cart_proj = wrf.get_cartopy(hgt)

    fig = plt.figure(figsize=(12,9),dpi=300)
    axe = plt.axes(projection=cart_proj)
    axe.set_title(figtle,loc='left',fontsize=title_font,fontdict=font) 
    axe.set_title('m',loc='right',fontsize=title_font,fontdict=font)
    
    axe.add_feature(cfeat.NaturalEarthFeature('physical', 'land', '50m',
        edgecolor='k',facecolor='none',linewidth=1.0))
    coast_shp2 = Reader(os.getenv("SHP_LIB")+"/cnmap/cnhimap.dbf").geometries()
    coastline2 = cfeat.ShapelyFeature(coast_shp2, ccrs.PlateCarree(), 
        edgecolor="k",facecolor="none")
    axe.add_feature(coastline2, linewidth=1.5,zorder=2)
    
    cont = axe.contourf(lons, lats, hgt, cnlevels, 
         transform=ccrs.PlateCarree(),cmap=ncmap,
         extend='both',norm=norm)
    plt.colorbar(cont, ax=axe, shrink=.8, pad=0.01)
    axe.set_xlim(wrf.cartopy_xlim(hgt))
    axe.set_ylim(wrf.cartopy_ylim(hgt))
    del lons, lats, hgt
    gc.collect()
    
    ds = xr.open_dataset(filname[1])
    lat2 = ds['lat_rho'].data
    lon2 = ds['lon_rho'].data
    axe.plot(np.hstack((lon2[0,:],lon2[:,-1],lon2[-1,::-1],lon2[:,0])),
        np.hstack((lat2[0,:],lat2[:,-1],lat2[-1,:],lat2[::-1,0])), linewidth=2.0, 
        linestyle='dashed', color='r', transform=ccrs.PlateCarree())

    gl = axe.gridlines(crs=ccrs.PlateCarree(),draw_labels=True, linewidth=1., 
        color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    #gl.xlocator = mticker.FixedLocator([-180, -45, 0, 45, 180])
    plt.savefig(figname,bbox_inches="tight",pad_inches=0.01)

def draw_wrf_var(nf,filname,varname,unit,cnlevels,ncmap):
    norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=ncmap.N, extend='both')
    
    d01 = Dataset(filname)
    time = datetime.strptime(np.datetime_as_string(
        wrf.getvar(d01,'Times').data),'%Y-%m-%dT%H:00:00.000000000')
    figname = '%s/%s_%s_%s.png'%(figdir,varname,nf,time.strftime('%Y%m%d%H'))
    if os.path.exists(figname):
        print('%s exists'%figname)
        return

    var = wrf.getvar(d01,varname)
    var.data = var.data-273.15
    print('%s : min = %f ; max = %f'%(varname,
        np.nanmin(var.data),np.nanmax(var.data)))
    lats, lons = wrf.latlon_coords(var)
    cart_proj = wrf.get_cartopy(var)

    fig = plt.figure(figsize=(8,6),dpi=300)
    axe = plt.axes(projection=cart_proj)
    axe.set_title('%s %s @ %s'%(nf,varname,time.strftime('%b %d %H:00')),
        loc='left',fontsize=title_font,fontdict=font) 
    axe.set_title(unit,loc='right',fontsize=title_font,fontdict=font)
    
    axe.add_feature(cfeat.NaturalEarthFeature('physical', 'land', '50m',
        edgecolor='k',facecolor='grey',linewidth=1.0))
    #coast_shp2 = Reader(os.getenv("SHP_LIB")+"/cnmap/cnhimap.dbf").geometries()
    #coastline2 = cfeat.ShapelyFeature(coast_shp2, ccrs.PlateCarree(), 
    #    edgecolor="k",facecolor="none")
    #axe.add_feature(coastline2, linewidth=1.5,zorder=2)
    
    cont = axe.contourf(lons, lats, var, cnlevels, 
         transform=ccrs.PlateCarree(),cmap=ncmap,
         extend='both',norm=norm)
    plt.colorbar(cont, ax=axe, shrink=.8, pad=0.01)
    axe.set_xlim(wrf.cartopy_xlim(var))
    axe.set_ylim(wrf.cartopy_ylim(var))
    del lons, lats, var
    gc.collect()
    
    ds = xr.open_dataset('/home/metctm1/array_hq133/data/s2s/wrfroms/'+
        'realtime/CP_free/njord_his_d01.20220725.nc')
    lat2 = ds['lat_rho'].data
    lon2 = ds['lon_rho'].data
    axe.plot(np.hstack((lon2[0,:],lon2[:,-1],lon2[-1,::-1],lon2[:,0])),
        np.hstack((lat2[0,:],lat2[:,-1],lat2[-1,:],lat2[::-1,0])), linewidth=2.0, 
        linestyle='dashed', color='k', transform=ccrs.PlateCarree())

    gl = axe.gridlines(crs=ccrs.PlateCarree(),draw_labels=True, linewidth=1., 
        color='k', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    plt.savefig(figname,bbox_inches="tight",pad_inches=0.01)
    del fig, axe, gl, ds
    gc.collect()

def animat(inputfile, gif_name): 
    fn_stream = subprocess.check_output("ls -rt %s"%inputfile, 
            shell=True).decode("utf-8")
    fn_list   = fn_stream.split()
    print(fn_list[0])
    print("filenumber : "+str(len(fn_list)))

    frames = []
    for itm in fn_list:
        frame = Image.open(itm)
        frames.append(frame)

    frames[0].save(gif_name, save_all=True, append_images=frames[1:],\
                duration = 250, loop=0, disposal=0)
    # duration: The time to display the current frame of the GIF, in milliseconds.
    # loop: The number of times the GIF should loop. 0 means that it will loop forever.
    subprocess.run("rm -f %s"%inputfile,shell=True)

if __name__=='__main__':
    main_run()

