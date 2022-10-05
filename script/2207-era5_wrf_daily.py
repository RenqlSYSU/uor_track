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

title_font=18
label_font=14
plt.rcParams["font.weight"] = "bold"
font = {'family': 'sans-serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black',
        }

indir = ['/home/metctm1/array_hq133/data/drv_fld/era5',
         '/home/metctm1/array_hq133/data/s2s/wrfonly/clim/',
         '/home/metctm1/array_hq133/data/s2s/wrfroms/clim/aegir_']
outdir = '/home/lzhenn/cooperate/data/cwrf_s2s'
figdir = '/home/lzhenn/cooperate/fig'
case = ['ERA5','WRF','WRFROMS']
years = range(2011,2021) #[2022,]#
lats = 15
latn = 50
lonl = 75
lonr = 137
ilat = np.arange(lats, latn+0.1, 0.25)
ilon = np.arange(lonl, lonr+0.1, 0.25)
ftime  = pd.date_range(start='2021-07-01 00',end='2021-08-05 00',
    freq='1D',closed=None)
print('%d days'%len(ftime))

def main_run():
    read_era5_grib('t2m','T2','K','sl')
    #for i in range(1,3):
    #    calc_wrf_interp('T2','K',3,indir[i],case[i])
    #draw_clim_tcc('T2','T2m',np.arange(16,33,1),np.arange(-2.4,2.41,0.3))
    #test_regionmask('T2','T2m')

def test_regionmask(varname,drawvar):
    var = [] 
    for ca in range(len(case)):
        if len(years)==1:
            ds = xr.open_dataset('%s/%s-%s_daily_%d.nc'%(outdir,case[ca],varname,years[0]))
        else:
            ds = xr.open_dataset('%s/%s-%s_daily_clim.nc'%(outdir,case[ca],varname))
        var.append(ds[varname].sel(time=ftime,lat=ilat,lon=ilon))
    # choose china grid
    #shf = '/disk/r074/lzhenn/tracacode/2205-PPOL-COOP/data/gadm36_CHN_0.shp'
    shf = '/home/metctm1/array/tracacode/UTILITY-2016/shp/cnmap/gadm36_CHN_1.shp'
    mask = regionmask.mask_geopandas(gpd.read_file(shf),var[0])
    var = [x.where(mask==5).mean(('lat','lon'))-273.15 for x in var]

    fig = plt.figure(figsize=(12,6),dpi=150)
    axe = fig.subplots(1,1)
    axe.set_title('Guangdong %s 2016'%drawvar,fontsize=title_font,fontdict=font)
    for nr in range(len(var)):
        axe.plot(ftime, var[nr], linewidth=2)
    
    axe.set_xlabel('',fontsize=label_font,fontdict=font)
    axe.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    plt.setp(axe.get_xticklabels(), rotation=15, ha="right")
    
    axe.set_ylabel('',fontsize=label_font,fontdict=font)
    
    axe.grid(True, which="both", color='grey', linestyle='--', linewidth=1)
    axe.legend(case)
    
    plt.savefig('%s/%s_ts.png'%(figdir,varname), 
        bbox_inches='tight',pad_inches=0.01)
    
def draw_clim_tcc(varname,drawvar,cnlev1,cnlev2):
    var = [] 
    for ca in range(len(case)):
        ds = xr.open_dataset('%s/%s-%s_daily_clim.nc'%(outdir,case[ca],varname))
        var.append(ds[varname].sel(time=ftime,lat=ilat,lon=ilon).data-273.15)
'''
    draw_3plot([var[1]-var[0],var[2]-var[0],var[2]-var[1]], 
        ['diff WRF-ERA5','diff WRFROMS-ERA5','diff WRFROMS-WRF'], 
        cnlev2, cmaps.BlueDarkRed18, '%s/cwrf_%s_diff.png'%(figdir,varname))
    draw_3plot([var[1]-var[0],var[2]-var[0],var[2]-var[1]], 
        ['diff WRF-ERA5','diff WRFROMS-ERA5','diff WRFROMS-WRF'], 
        np.arange(-0.8,0.81,0.1), cmaps.BlueDarkRed18, '%s/cwrf_%s_diff2.png'%(figdir,varname))
    draw_3plot(var, ['%s %s'%(drawvar,nc) for nc in case],cnlev1, 
        cmaps.precip2_17lev, '%s/cwrf_%s_clim.png'%(figdir,varname))
'''

def draw_3plot(var,figtitle,cnlevels,fcolors,figfile):
    lat_sp = 10
    nrow = len(var) #2 #
    ncol = 1 #6 #
    norm = colors.BoundaryNorm(boundaries=cnlevels, 
        ncolors=fcolors.N,extend='both')
    
    fig = plt.figure(figsize=(12,12),dpi=300)
    ax = fig.subplots(nrow, ncol)
    for nr in range(0,len(var),1):
        print(var[nr].min())
        print(var[nr].max())
        axe = ax[nr]
        axe.set_title(figtitle[nr],fontsize=title_font)
        cont = axe.contourf(ftime, ilat, var[nr].transpose(), cnlevels, 
             cmap=fcolors,extend='both',norm=norm)
        plt.colorbar(cont, ax=axe, shrink=.9, pad=0.01)
        
        axe.set_xlabel('',fontsize=label_font,fontdict=font)
        axe.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
        plt.setp(axe.get_xticklabels(), rotation=15, ha="right")
        
        axe.set_ylabel('',fontsize=label_font,fontdict=font)
        axe.set_yticks(np.arange(lats,latn,lat_sp))
        axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))

    plt.savefig(figfile, bbox_inches='tight',pad_inches=0.01)
    
def read_era5_grib(varname,new_var,unit,suffix):
    if len(years)==1:
        outfile = '%s/ERA5-%s_daily_%d.nc'%(outdir,new_var,years[0])
    else:
        outfile = '%s/ERA5-%s_daily_clim.nc'%(outdir,new_var)
    if os.path.exists(outfile):
        print('%s exists'%outfile)
        return
    else:
        print('handle %s'%outfile)
    
    for year in years:
        fn_stream = subprocess.check_output(
            'ls %s/%d0[78]*-%s.grib'%(
            indir[0],year,suffix), shell=True).decode('utf-8')
        fn_list = fn_stream.split()[0:len(ftime)]
        term = []
        for itm in range(len(fn_list)):
            print(fn_list[itm])
            ds = xr.open_dataset(fn_list[itm],engine='cfgrib',backend_kwargs={
                'filter_by_keys': {'typeOfLevel': 'surface'},'indexpath':''})
            term.append(ds[varname].mean('time').data.tolist())
            #term.append(ds[varname].sel(
            #    longitude=ilon).mean(('time','longitude')).data.tolist())
        if year == years[0]:
            var = np.array(term)
            print('var shape %s'%(str(var.shape)))
        else:
            var = var + np.array(term)
    print('var shape %s'%(str(var.shape)))
    da = xr.DataArray(np.divide(var.data,len(years)),coords=[ftime,
        ds.latitude,ds.longitude], dims=['time','lat','lon'])
    da.attrs['units']=unit
    ds = da.to_dataset(name=new_var)
    ds.to_netcdf(outfile,"w")

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
            'ls %s%d070100/wrfout_d01_%d-07-*'%(
            path,year,year), shell=True).decode('utf-8')
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
