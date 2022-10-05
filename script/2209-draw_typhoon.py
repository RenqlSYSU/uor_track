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
paths = ['/home/lzhenn/cooperate/data/case_study/coupled/2018091200',
         '/home/lzhenn/cooperate/data/case_study/coupled/2018091200pgw_nudging',
         '/home/lzhenn/cooperate/data/case_study/coupled/2018091200pgw_free']
case = ['CTRL','2090_PGW_NG','2090_PGW_FREE']
clrs = ['blue','red','green']
#ftime  = pd.date_range(start='2018-09-15 06',end='2018-09-17 00',
#    freq='1H',closed=None)
ftime  = pd.date_range(start='2018-09-16 00',end='2018-09-16 07',
    freq='3H',closed=None)
#ftime  = pd.date_range(start='2018-09-14 00',end='2018-09-17 01',
#    freq='3D',closed=None)
domain = 'd02'
print(ftime)

def main_run():
    #for nt in ftime:
    #    filname = '%s/wrfout_%s_%s'%(paths[2],domain,nt.strftime("%Y-%m-%d_%H:00:00"))
    #    draw_wrf_var('2090_PGW_FREE-CTRL',filname,['wspd10',],'m/s',
    #        np.arange(-8,8.1,1), cmaps.BlueDarkRed18)
    #    draw_wrf_var('2090_PGW_FREE-CTRL',filname,['slp',],'hPa',
    #        np.arange(-8,8.1,1), cmaps.BlueDarkRed18)

    #cnlevels = np.arange(20,31.6,0.1)
    #ncmap = colors.ListedColormap(cmaps.rainbow(
    #    np.linspace(0,1,len(cnlevels)+1))) 
    #cnlevels = [0.1, 0.5, 1, 3, 5, 10,15,20, 30, 40, 60, 80, 100, 120, 150, 200, 250] #24h accumulated preci
    #cnlevels = [0.1, 1, 3, 5, 10, 15, 20, 30, 50, 70, 100, 150, 200, 300, 400, 500, 600] #24h accumulated preci
    for nc in range(0,len(case)):
        filname = '%s/wrfout_%s_%s'%(paths[nc],domain,ftime[-1].strftime("%Y-%m-%d_%H:00:00"))
        draw_wrf_var(case[nc],filname,['RAINNC',],'mm',cnlevels,cmaps.precip2_17lev)
        #calc_minslp(case[nc],paths[nc])
    #draw_track()
    #draw_ts(4,'Minimum Sea Level Pressure (hPa)','%s/mslp_%s.png'%(figdir,domain))
    #draw_ts(5,'Maximum Wind Speed (m/s)','%s/wspd_%s.png'%(figdir,domain))

def draw_ts(ind,title,figname):
    fig = plt.figure(figsize=(12,4),dpi=150) #width,height
    axe = fig.subplots(1,1)
    axe.set_title(title,fontsize=title_font,fontdict=font)
    
    for nc in range(len(case)):
        data = np.loadtxt('/home/lzhenn/cooperate/data/%s_%s_minslp.txt'%(
            case[nc],domain),usecols=(0,ind))
        axe.plot(pd.to_datetime(['%d'%x for x in data[:,0]],format='%Y%m%d%H'),
            data[:,1],'o-',linewidth=1.0,color=clrs[nc],markersize=2.5,
            label=case[nc])
    data = np.loadtxt('/home/metctm1/array_hq86/data/1911-COAWST/hko.trck.mangkhut',
        skiprows=1,usecols=(0,ind))
    axe.plot(pd.to_datetime(['%d'%x for x in data[:,0]],format='%Y%m%d%H'),
        data[:,1],'o-',linewidth=1.0,color='black',markersize=2.5,
        label='HKO bestTrack')
    data = np.loadtxt('/home/metctm1/array_hq86/data/1911-COAWST/cma.trck.mangkhut',
        skiprows=1,usecols=(0,ind))
    axe.plot(pd.to_datetime(['%d'%x for x in data[:,0]],format='%Y%m%d%H'),
        data[:,1],'o-',linewidth=1.0,color='grey',markersize=2.5,
        label='CMA bestTrack')
    plt.legend(loc='best', fontsize=label_font)
    
    axe.set_xlim(ftime[0],ftime[-1])
    if ind==4 :
        axe.set_ylim(940,1000)
    axe.set_ylabel('',fontsize=label_font,fontdict=font)
    axe.set_xlabel('time',fontsize=label_font,fontdict=font)
    axe.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d %H'))
    plt.setp(axe.get_xticklabels(), rotation=15, ha="right")
    plt.savefig(figname,bbox_inches="tight",pad_inches=0.01)

def draw_track():
    title = 'Tracks of Mangkhut(2018)'
    figname = '%s/Mangkhut_track.png'%(figdir)
    
    filname = '%s/wrfout_%s_%s'%(paths[1],domain,ftime[0].strftime("%Y-%m-%d_%H:00:00"))
    ncfile = Dataset(filname)
    var    = wrf.getvar(ncfile, 'slp')
    lats, lons = wrf.latlon_coords(var)
    cart_proj = wrf.get_cartopy(var)

    fig = plt.figure(figsize=(8,6),dpi=150)
    axe = plt.axes(projection=cart_proj)
    axe.set_title(title, fontsize=title_font, fontdict=font)
    
    MAP_RES = '50m'
    axe.coastlines(MAP_RES, linewidth=0.8)
    axe.add_feature(cfeat.OCEAN.with_scale(MAP_RES))
    axe.add_feature(cfeat.LAND.with_scale(MAP_RES))
    axe.add_feature(cfeat.LAKES.with_scale(MAP_RES))
    coast_shp2 = Reader(os.getenv("SHP_LIB")+"/cnmap/cnhimap.dbf").geometries()
    coastline2 = cfeat.ShapelyFeature(coast_shp2, ccrs.PlateCarree(), 
        edgecolor="k",facecolor="none")
    axe.add_feature(coastline2, linewidth=0.5,zorder=2)
    axe.set_xlim(wrf.cartopy_xlim(var))
    axe.set_ylim(wrf.cartopy_ylim(var))
    del lons, lats, var, ncfile
    gc.collect()

    if domain == 'd01':
        filname = '%s/wrfout_d02_%s'%(paths[1],ftime[0].strftime("%Y-%m-%d_%H:00:00"))
        ncfile = Dataset(filname)
        var    = wrf.getvar(ncfile, 'slp')
        lat2, lon2 = wrf.latlon_coords(var)
        axe.plot(np.hstack((lon2[0,:],lon2[:,-1],lon2[-1,::-1],lon2[:,0])),
            np.hstack((lat2[0,:],lat2[:,-1],lat2[-1,:],lat2[::-1,0])), linewidth=1.5, 
            linestyle='dashed', color='r', transform=ccrs.PlateCarree())
    
    data = np.loadtxt('/home/metctm1/array_hq86/data/1911-COAWST/hko.trck.mangkhut',
        skiprows=1,usecols=(2,3))
    axe.plot(data[:,1]/10.0,data[:,0]/10.0,'o-',linewidth=1.0,color='black',markersize=2.5,
        transform=ccrs.PlateCarree(),label='HKO bestTrack')
    data = np.loadtxt('/home/metctm1/array_hq86/data/1911-COAWST/cma.trck.mangkhut',
        skiprows=1,usecols=(2,3))
    axe.plot(data[:,1]/10.0,data[:,0]/10.0,'o-',linewidth=1.0,color='grey',markersize=2.5,
        transform=ccrs.PlateCarree(),label='CMA bestTrack')
    for nc in range(len(case)):
        data = np.loadtxt('/home/lzhenn/cooperate/data/%s_%s_minslp.txt'%(
            case[nc],domain),usecols=(2,3))
        axe.plot(data[:,1],data[:,0],'o-',linewidth=1.0,color=clrs[nc],markersize=2.5,
            transform=ccrs.PlateCarree(),label=case[nc])
    plt.legend(loc='best', fontsize=label_font)
    
    gl = axe.gridlines(crs=ccrs.PlateCarree(),draw_labels=True, linewidth=1., 
        color='k', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    plt.savefig(figname,bbox_inches="tight",pad_inches=0.01)
    del fig, axe, gl
    gc.collect()
    
def calc_minslp(prefix, path):
    outfile = '/home/lzhenn/cooperate/data/%s_%s_minslp.txt'%(prefix,domain)
    if os.path.exists(outfile):
        print('%s exists'%outfile)
        return
    
    ff = open(outfile,'w')
    for time in ftime:
        filname = '%s/wrfout_%s_%s'%(path,domain,time.strftime("%Y-%m-%d_%H:00:00"))
        print('read %s'%filname)
        ncfile = Dataset(filname)
        var    = wrf.getvar(ncfile, 'slp')
        ind = np.unravel_index(np.argmin(var.data), var.shape)
        lat, lon = wrf.latlon_coords(var)
        wspd   = wrf.getvar(ncfile, 'wspd_wdir10')[0].data
        ff.write('%s %s %f %f %f %f\n'%(time.strftime("%Y%m%d%H"),domain,
            lat[ind],lon[ind],var[ind],np.max(wspd)))
    ff.close()

def draw_wrf_var(nf,filname,varname,unit,cnlevels,ncmap):
    norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=ncmap.N, extend='both')
    
    d01 = Dataset(filname)
    time = datetime.strptime(np.datetime_as_string(
        wrf.getvar(d01,'Times').data),'%Y-%m-%dT%H:00:00.000000000')
    figname = '%s/%s_%s_%s_%s.png'%(figdir,varname[0],domain,nf,time.strftime('%Y%m%d%H'))
    if os.path.exists(figname):
        print('%s exists'%figname)
        return
    
    if varname[0] == 'wspd10':
        var = wrf.getvar(d01,'wspd_wdir10')[0]
        print(var)
    else:
        var = wrf.getvar(d01,varname[0])

    if varname[0] in ['SST',]:
        var.data = var.data-273.15
    if varname[0] in ['RAINNC',]:
        d011 = Dataset('%s-%s'%(filname.split('-',1)[0],
            ftime[0].strftime("%m-%d_%H:00:00")))
        var.data = var.data + wrf.getvar(d01,'RAINC').data - wrf.getvar(
            d011,'RAINNC').data - wrf.getvar(d011,'RAINC').data
        varname[0] = 'preci'
        print('%s %s: %f, %f'%(nf,varname[0], 
            np.ma.array(var.data,mask=(var.data<0.1)).mean(),
            var.data.mean()))
        del d011
    if nf == ['2090_PGW_FREE-CTRL','2090_PGW-CTRL']:
        d011 = Dataset('%s/wrfout_%s_%s'%(paths[0],domain,
            time.strftime("%Y-%m-%d_%H:00:00")))
        if varname[0] == 'wspd10':
            var.data = var.data - wrf.getvar(d011,'wspd_wdir10')[0].data
        else:
            var.data = var.data - wrf.getvar(d011,varname[0]).data
        del d011
        
    print('%s : min = %f ; max = %f'%(varname[0],
        np.nanmin(var.data),np.nanmax(var.data)))
    lats, lons = wrf.latlon_coords(var)
    cart_proj = wrf.get_cartopy(var)

    fig = plt.figure(figsize=(8,6),dpi=150)
    axe = plt.axes(projection=cart_proj)
    axe.set_title('%s %s @ %s'%(nf,varname[0],time.strftime('%b %d %H:00')),
        loc='left',fontsize=title_font,fontdict=font) 
    axe.set_title(unit,loc='right',fontsize=title_font,fontdict=font)
    
    axe.add_feature(cfeat.NaturalEarthFeature('physical', 'land', '50m',
        edgecolor='k',facecolor='None',linewidth=1.0))
    coast_shp2 = Reader(os.getenv("SHP_LIB")+"/cnmap/cnhimap.dbf").geometries()
    coastline2 = cfeat.ShapelyFeature(coast_shp2, ccrs.PlateCarree(), 
        edgecolor="k",facecolor="none")
    axe.add_feature(coastline2, linewidth=1.5,zorder=2)
    
    cont = axe.contourf(lons, lats, var, cnlevels, 
         transform=ccrs.PlateCarree(),cmap=ncmap,
         extend='both',norm=norm)
    plt.colorbar(cont, ax=axe, shrink=.9, pad=0.01)
    axe.set_xlim(wrf.cartopy_xlim(var))
    axe.set_ylim(wrf.cartopy_ylim(var))
    del lons, lats, var
    gc.collect()
    
    gl = axe.gridlines(crs=ccrs.PlateCarree(),draw_labels=True, linewidth=1., 
        color='k', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    plt.savefig(figname,bbox_inches="tight",pad_inches=0.01)
    del fig, axe, gl
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

