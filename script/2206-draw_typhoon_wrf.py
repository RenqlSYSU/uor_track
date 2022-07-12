import xarray as xr
import pandas as pd
import numpy as np
import scipy.io as scio
import subprocess
import os 
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.colors import ListedColormap
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.io.shapereader import Reader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from PIL import Image
import cmaps
from netCDF4 import Dataset
from wrf import (to_np, getvar, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, interplevel)

title_font=12
label_font=10
plt.rcParams["font.weight"] = "bold"
font = {'family': 'sans-serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black',
        }

figdir = "/home/lzhenn/cooperate/fig"
wrfpath = '/home/lzhenn/array127/data/oriental/2021090900'
swnpath = '/home/lzhenn/cooperate/data/case_study/chanthu/2021090900'
swngird = '/home/lzhenn/array74/workspace/calypso_pipeline/domaindb/oriental/roms_d01.nc'
ftime  = pd.date_range(start='2021-09-10 00', end='2021-09-18 00' \
    , freq='1H', closed=None)
lat_sp = 2.0 #5.0
lon_sp = 2.0 #10.0

def main_run():
    #calc_minslp()
    ds = xr.open_dataset(swngird)
    lat = ds['lat_rho'].data
    lon = ds['lon_rho'].data
    i = 70
    ''' 
    for i in range(len(ftime)):
        print(ftime[i].strftime("%Y-%m-%d %HZ"))
        draw_map_swan_hsig(i,ftime[i],lat,lon)
        draw_map_swan_hswell(i,ftime[i],lat,lon)
        draw_map_wrf_slp(i,ftime[i])
    '''
    animat('%s/slp_*.png'%figdir,'%s/slp_uv10.gif'%figdir)
    animat('%s/hsig_*.png'%figdir,'%s/hsig_wave.gif'%figdir)
    animat('%s/hswell_*.png'%figdir,'%s/hswell.gif'%figdir)

def draw_map_swan_hswell(i,time,lat,lon):
    title = 'Hswell @ %s'%(ftime[i].strftime("%Y-%m-%d %HZ"))
    figname = '%s/hswell_%s.png'%(figdir,ftime[i].strftime("%Y-%m-%d_%H"))
    data = scio.loadmat('%s/hswell_%s.00.mat'%(swnpath,time.strftime("%Y%m%d")))
    var  = data['Hswell_%s0000'%time.strftime("%Y%m%d_%H")] 
    
    fig = plt.figure(figsize=(12,9),dpi=300)
    axe = plt.axes(projection=ccrs.PlateCarree())
    axe.set_title(title, fontsize=title_font, fontdict=font)
    axe.add_feature(cfeat.COASTLINE.with_scale('10m'), linewidth=1,color='k')
    coast_shp1 = Reader('/home/lzhenn/njord_pipeline/postprocess/shp/cnhimap.dbf').geometries()
    coastline1 = cfeat.ShapelyFeature(coast_shp1, ccrs.PlateCarree(), edgecolor="black", facecolor="none")
    axe.add_feature(coastline1, linewidth=0.8,zorder=2)
    
    cnlevels = np.arange(0.5,8.5,0.5) 
    ncmap = cmaps.precip3_16lev
    norm  = colors.BoundaryNorm(boundaries=cnlevels, ncolors=ncmap.N,extend="both")
    cont = axe.pcolormesh(lon, lat, var, zorder=0, 
             transform=ccrs.PlateCarree(),cmap=ncmap,norm=norm)
    plt.colorbar(cont, ax=axe, shrink=.7)

    draw_track(i,axe,lon,lat)
    plt.savefig(figname,bbox_inches="tight",pad_inches=0.)

def draw_map_swan_hsig(i,time,lat,lon):
    q_mis=20 # wind vector plotting every q_miss grid
    title = 'Hsig & wave @ %s'%(ftime[i].strftime("%Y-%m-%d %HZ"))
    figname = '%s/hsig_%s.png'%(figdir,ftime[i].strftime("%Y-%m-%d_%H"))
    data = scio.loadmat('%s/hsig_%s.00.mat'%(swnpath,time.strftime("%Y%m%d")))
    var  = data['Hsig_%s0000'%time.strftime("%Y%m%d_%H")] 
    
    data = scio.loadmat('%s/lwavp_%s.00.mat'%(swnpath,time.strftime("%Y%m%d")))
    lwavp= data['Lwavp_%s0000'%time.strftime("%Y%m%d_%H")] 
    data = scio.loadmat('%s/tps_%s.00.mat'%(swnpath,time.strftime("%Y%m%d")))
    tps  = data['TPsmoo_%s0000'%time.strftime("%Y%m%d_%H")] 
    data = scio.loadmat('%s/pdir_%s.00.mat'%(swnpath,time.strftime("%Y%m%d")))
    pdir = data['PkDir_%s0000'%time.strftime("%Y%m%d_%H")]+180.0 
    uvar = np.sin(pdir*np.pi/180.0) * lwavp / tps 
    vvar = np.cos(pdir*np.pi/180.0) * lwavp / tps 
    
    fig = plt.figure(figsize=(12,9),dpi=300)
    axe = plt.axes(projection=ccrs.PlateCarree())
    axe.set_title(title, fontsize=title_font, fontdict=font)
    axe.add_feature(cfeat.COASTLINE.with_scale('10m'), linewidth=1,color='k')
    coast_shp1 = Reader('/home/lzhenn/njord_pipeline/postprocess/shp/cnhimap.dbf').geometries()
    coastline1 = cfeat.ShapelyFeature(coast_shp1, ccrs.PlateCarree(), edgecolor="black", facecolor="none")
    axe.add_feature(coastline1, linewidth=0.8,zorder=2)
    
    cnlevels = np.arange(0.5,8.5,0.5) 
    ncmap = cmaps.precip3_16lev
    norm  = colors.BoundaryNorm(boundaries=cnlevels, ncolors=ncmap.N,extend="both")
    cont = axe.pcolormesh(lon, lat, var, zorder=0, 
             transform=ccrs.PlateCarree(),cmap=ncmap,norm=norm)
    plt.colorbar(cont, ax=axe, shrink=.7)

    quv = axe.quiver(lon[::q_mis,::q_mis],lat[::q_mis,::q_mis],
           uvar[::q_mis,::q_mis],vvar[::q_mis,::q_mis],zorder=2,
           pivot='mid',units='inches',scale=60,scale_units='inches',color="dimgray",
           width=0.02,headwidth=3,headlength=4.5,transform=ccrs.PlateCarree())
    plt.quiverkey(quv, 1.1, 0.1, 20, r'$20 m/s$', labelpos='N',
           coordinates='axes')
    
    draw_track(i,axe,lon,lat)
    plt.savefig(figname,bbox_inches="tight",pad_inches=0.)

def draw_map_wrf_slp(i,time):
    q_mis=12 # wind vector plotting every q_miss grid
    title = 'SLP(hPa) & UV10 @ %s'%(time.strftime("%Y-%m-%d %HZ"))
    figname = '%s/slp_%s.png'%(figdir,time.strftime("%Y-%m-%d_%H"))
    ncfile = Dataset('%s/wrfout_d02_%s'%(wrfpath,time.strftime("%Y-%m-%d_%H:00:00")))
    var    = getvar(ncfile, 'slp')
    uvar   = getvar(ncfile, 'U10').data
    vvar   = getvar(ncfile, 'V10').data
    lat, lon = latlon_coords(var)
    cart_proj = get_cartopy(var)
    
    fig = plt.figure(figsize=(12,9),dpi=300)
    axe = plt.axes(projection=cart_proj)
    axe.set_title(title, fontsize=title_font, fontdict=font)
    axe.add_feature(cfeat.COASTLINE.with_scale('10m'), linewidth=1,color='k')
    coast_shp1 = Reader('/home/lzhenn/njord_pipeline/postprocess/shp/cnhimap.dbf').geometries()
    coastline1 = cfeat.ShapelyFeature(coast_shp1, ccrs.PlateCarree(), edgecolor="black", facecolor="none")
    axe.add_feature(coastline1, linewidth=0.8,zorder=2)
    
    cnlevels = np.arange(975,1023,3) 
    ncmap = ListedColormap(cmaps.precip3_16lev(range(0,17,1))[::-1]) 
    norm  = colors.BoundaryNorm(boundaries=cnlevels, ncolors=ncmap.N,extend="both")
    cont = axe.pcolormesh(lon, lat, var, zorder=0, 
             transform=ccrs.PlateCarree(),cmap=ncmap,norm=norm)
    plt.colorbar(cont, ax=axe, shrink=.7)

    quv = axe.quiver(to_np(lon[::q_mis,::q_mis]),to_np(lat[::q_mis,::q_mis]),
           uvar[::q_mis,::q_mis],vvar[::q_mis,::q_mis],zorder=2,
           pivot='mid',units='inches',scale=30,scale_units='inches',color="dimgray",
           width=0.02,headwidth=3,headlength=4.5,transform=ccrs.PlateCarree())
    plt.quiverkey(quv, 1.1, 0.1, 10, r'$10 m/s$', labelpos='N',
           coordinates='axes')
    
    draw_track(i,axe,lon,lat)
    axe.set_xlim(cartopy_xlim(var))
    axe.set_ylim(cartopy_ylim(var))
    plt.savefig(figname,bbox_inches="tight",pad_inches=0.)

def calc_minslp():
    ff = open('/home/lzhenn/cooperate/data/chanthu_wrf_minslp.txt','w')
    for time in ftime:
        filname = '%s/wrfout_d02_%s'%(wrfpath,time.strftime("%Y-%m-%d_%H:00:00"))
        ncfile = Dataset(filname)
        var    = getvar(ncfile, 'slp')
        ind = np.unravel_index(np.argmin(var.data), var.shape)
        lat, lon = latlon_coords(var)
        if var[ind] < 1000:
            ff.write('%s %f %f %f \n'%(time.strftime("%Y%m%d%H"),lat[ind],lon[ind],var[ind]))
        else:
            ff.write('%s %f %f %f \n'%(time.strftime("%Y%m%d%H"),0.0,180.0,var[ind]))
    ff.close()

def draw_track(i,axe,lon,lat):
    data = np.loadtxt('/home/lzhenn/cooperate/data/chanthu_bst_track.txt',
        skiprows=1,usecols=(2,3))
    tlat = data[:,0]/10.0
    tlon = data[:,1]/10.0
    axe.plot(tlon,tlat,'o-',linewidth=2.5, color='black', 
        transform=ccrs.PlateCarree())
    
    data = np.loadtxt('/home/lzhenn/cooperate/data/chanthu_wrf_minslp.txt',
        usecols=(1,2))
    tlat = data[0:i+1,0]
    tlon = data[0:i+1,1]
    axe.plot(tlon,tlat,'-',linewidth=2.5, color='m', 
        transform=ccrs.PlateCarree())
    
    axe.set_xlim(lon[0,0],lon[-1,-1])
    axe.set_ylim(lat[0,0],lat[-1,-1])
    axe.set_xticks(np.arange(np.ceil(lon[0,0]),np.floor(lon[0,-1]),lon_sp), crs=ccrs.PlateCarree())
    axe.set_yticks(np.arange(np.ceil(lat[0,0]),np.floor(lat[-1,0]),lat_sp), crs=ccrs.PlateCarree())
    axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
    axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))

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

