import xarray as xr
import pandas as pd
import numpy as np
from scipy import stats
import gc #garbage collector
from scipy.interpolate import griddata
import os, subprocess, wrf
from multiprocessing import Pool
from netCDF4 import Dataset
import wrf
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.dates as mdates
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.io.shapereader import Reader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmaps, regionmask
import geopandas as gpd

title_font=14
label_font=10
plt.rcParams["font.weight"] = "bold"
font = {'family': 'sans-serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black',
        }

paths = ['/home/metctm1/array_hq133/data/s2s/wrfonly/realtime',
         '/home/metctm1/array_hq133/data/s2s/wrfroms/realtime']
indir = ['%s/WRF_NG'%paths[0],
         '%s/WRF_free'%paths[0],
         '%s/CFS_WRF_free'%paths[0],
         '%s/CP_NG'%paths[1],
         '%s/CP_free'%paths[1],
         '%s/CFS_CP_free'%paths[1],
         '%s/CFSA_HYCOM'%paths[1],
         '%s/ERA5_CFSO'%paths[1],
         '/home/metctm1/array_hq133/data/drv_fld/era5']
outdir = '/home/lzhenn/cooperate/data/cwrf_s2s'
figdir = '/home/lzhenn/cooperate/fig'
case = [x.split('/')[-1] for x in indir] 
case[-1]='ERA5'
years = [2022,] #range(2011,2021)
lats = 15
latn = 50
lonl = 75
lonr = 137
ilat = np.arange(lats, latn+0.1, 0.25)
ilon = np.arange(lonl, lonr+0.1, 0.25)
ftime  = pd.date_range(start='2021-07-01 00',end='2021-08-04 23',
    freq='1D',closed=None)
print('%d days'%len(ftime))

def main_run():
    #for i in range(1,len(indir)):
    #    calc_wrf_interp('T2','K',3,indir[i],case[i])
    draw_3x3plot('T2',np.arange(-2.4,2.41,0.3),cmaps.BlueDarkRed18)
    #regionmask_time_series('T2','T2m')

def draw_3x3plot(varname,cnlevels,fcolors):
    #ds = xr.open_dataset('%s/ERA5-%s_daily_2022.nc'%(outdir,varname))
    ds = xr.open_dataset('%s/ERA5-%s_daily_clim.nc'%(outdir,varname))
    mean = ds[varname].sel(time=ftime,lat=ilat,lon=ilon).mean('time').data
    var = [] 
    for ca in range(len(case)):
        ds = xr.open_dataset('%s/%s-%s_daily_%d.nc'%(outdir,case[ca],varname,years[0]))
        var.append(ds[varname].sel(time=ftime,lat=ilat,lon=ilon).mean('time').data-mean)
    
    lat_sp = 10
    lon_sp = 20 #60 #
    nrow = 3 #
    ncol = 3 #6 
    bmlo = 0.35
    prefix = [chr(i) for i in range(97,110)]
    norm = colors.BoundaryNorm(boundaries=cnlevels, 
        ncolors=fcolors.N,extend='both')
    
    fig = plt.figure(figsize=(12,12),dpi=300)
    ax = fig.subplots(nrow, ncol, subplot_kw=dict(
        projection=ccrs.PlateCarree(central_longitude=180.0)))
    for nr in range(0,len(var),1):
        print('min %f ; max %f'%(np.nanmin(var[nr]),np.nanmax(var[nr])))
        a,b = divmod(nr,3)
        axe = ax[a][b]
        axe.set_title('(%s) %s'%(prefix[nr],case[nr]),fontsize=title_font,fontdict=font)
        axe.add_feature(cfeat.NaturalEarthFeature('physical', 'land', '50m',
            edgecolor='k',facecolor='none',linewidth=0.8))
        coast_shp2 = Reader(os.getenv("SHP_LIB")+"/cnmap/cnhimap.dbf").geometries()
        coastline2 = cfeat.ShapelyFeature(coast_shp2, ccrs.PlateCarree(), 
            edgecolor="k",facecolor="none")
        axe.add_feature(coastline2, linewidth=0.8,zorder=2)

        cont = axe.contourf(ilon, ilat, var[nr], cnlevels, 
             transform=ccrs.PlateCarree(),cmap=fcolors,
             extend='both',norm=norm)
       
        if b==0:
            axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
            axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
        if a==2:
            axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
            axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))
    
    position = fig.add_axes([0.2, bmlo+0.005, 0.7, 0.01]) #left, bottom, width, height
    cb = plt.colorbar(cont, cax=position ,orientation='horizontal')#, shrink=.9)
    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig('%s/%s_diff.png'%(figdir,varname), 
        bbox_inches='tight',pad_inches=0.01)
    
def regionmask_time_series(varname,drawvar):
    lwidth = 9*[2] + [4]
    lcolor = 2*['r','b','g'] + ['c','m','k','k'] 
    lpattn = 3*['--'] + 7*['-']

    var = [] 
    for ca in range(len(case)):
        ds = xr.open_dataset('%s/%s-%s_daily_%d.nc'%(outdir,case[ca],varname,years[0]))
        var.append(ds[varname].sel(time=ftime,lat=ilat,lon=ilon))
    ds = xr.open_dataset('%s/ERA5-%s_daily_clim.nc'%(outdir,varname))
    var.append(ds[varname].sel(time=ftime,lat=ilat,lon=ilon))
    case.append('ERA5_clim')
    
    # choose china grid
    shf = '/disk/r074/lzhenn/tracacode/2205-PPOL-COOP/data/gadm36_CHN_0.shp'
    #shf = '/home/metctm1/array/tracacode/UTILITY-2016/shp/cnmap/gadm36_CHN_1.shp'
    mask = regionmask.mask_geopandas(gpd.read_file(shf),var[0])
    var = [x.where(mask==0).mean(('lat','lon'))-273.15 for x in var]

    fig = plt.figure(figsize=(12,6),dpi=150)
    axe = fig.subplots(1,1)
    axe.set_title('China %s 2022'%drawvar,fontsize=title_font,fontdict=font)
    for nr in range(len(var)):
        axe.plot(ftime, var[nr], linewidth=lwidth[nr], color=lcolor[nr],
            linestyle=lpattn[nr])
    
    axe.set_xlabel('',fontsize=label_font,fontdict=font)
    axe.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    plt.setp(axe.get_xticklabels(), rotation=15, ha="right")
    
    axe.set_ylabel('',fontsize=label_font,fontdict=font)
    
    axe.grid(True, which="both", color='grey', linestyle='--', linewidth=1)
    axe.legend(case)
    
    plt.savefig('%s/%s_ts.png'%(figdir,varname), 
        bbox_inches='tight',pad_inches=0.01)
    
def calc_wrf_interp(varname,unit,dh,path,casename):
    if len(years)==1:
        outfile = '%s/%s-%s_daily_%d.nc'%(outdir,casename,varname,years[0])
    else:
        outfile = '%s/%s-%s_daily_clim.nc'%(outdir,casename,varname)
    if os.path.exists(outfile):
        print('%s exists'%outfile)
        return
    else:
        print('handle %s'%outfile)
    
    for year in years:
        fn_stream = subprocess.check_output(
            'ls %s/wrfout_d01_%d-0[78]-*'%(
            path,year), shell=True).decode('utf-8')
        fn_list = fn_stream.split()
        wrf_list=[Dataset(itm) for itm in fn_list[::dh]]
        print('%d total file %d; opened file %d'%(year, 
            len(fn_list),len(wrf_list)))
        
        z = wrf.getvar(wrf_list, varname, timeidx=wrf.ALL_TIMES, 
            method='cat')
        if year == years[0]:
            term = z.groupby(z.Time.dt.dayofyear).mean('Time')
            xlat = wrf.getvar(wrf_list[0],'XLAT')
            xlon = wrf.getvar(wrf_list[0],'XLONG')
            print(term)
        else:
            term.data = term.data + z.groupby(z.Time.dt.dayofyear
                ).mean('Time').data 
        del wrf_list
        gc.collect()

    term.data = term.data/len(years)
    var = xr.apply_ufunc(grid_interp,term,xlat,xlon,
        input_core_dims=[['south_north','west_east'],
        ['south_north','west_east'],['south_north','west_east']],
        output_core_dims=[['lat','lon']],
        vectorize=True, dask="parallelized",output_dtypes=['float64'],)
    print(var)
                        
    da = xr.DataArray(var.data,
        coords=[ftime,ilat,ilon],dims=['time','lat','lon'])
    da.attrs['units'] = unit
    ds = da.to_dataset(name=varname)
    ds.to_netcdf(outfile,"w")

def grid_interp(var,lat,lon):
    '''
    All input variables are 2D numpy array
    1. compress all the variables into one dimension
    2. use scipy.interpolate.griddata to interp 
    '''
    dst_lon, dst_lat = np.meshgrid(ilon,ilat)
    points = np.column_stack((lat.flatten(),lon.flatten())) #(nlat*nlon,2)
    dst_points = np.column_stack((dst_lat.flatten(),
                                  dst_lon.flatten())) #(nlat*nlon,2)
    dst_var = griddata(points, var.flatten(), dst_points,
        method='linear')
    return dst_var.reshape(dst_lat.shape)

if __name__=='__main__':
    main_run()
