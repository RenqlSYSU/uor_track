import xarray as xr
import pandas as pd
import numpy as np
import gc #garbage collector
from scipy.interpolate import griddata
import os, subprocess, wrf
from multiprocessing import Pool
from netCDF4 import Dataset
import wrf
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.io.shapereader import Reader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmaps

title_font=18
label_font=14
plt.rcParams["font.weight"] = "bold"
font = {'family': 'sans-serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black',
        }

indir = ['/home/metctm1/cmip-wrfltm-arch/era5',
         '/home/metctm1/array_hq133/data/s2s/wrfonly/clim/',
         '/home/metctm1/array_hq133/data/s2s/wrfroms/clim/aegir_']
outdir = '/home/lzhenn/cooperate/data/cwrf_s2s'
figdir = '/home/lzhenn/cooperate/fig'
case = ['ERA5','WRF','WRFROMS']
lats = 15
latn = 50
lonl = 75
lonr = 137
ilat = np.arange(lats, latn+0.1, 0.25)
ilon = np.arange(lonl, lonr+0.1, 0.25)

def main_run():
#    for i in range(1,3):
        for year in range(2011,2021):
            read_era5_grib(year,'t2m','K','sl')
#            calc_wrf_interp(year,'T2',1,indir[i],case[i])

def draw_clim_tcc():
    ds = xr.open_mfdataset('%s/%s_%s_20*.nc'%(outdir,casename,varname),
       combine='nested', concat_dim='time', parallel=True)
    var = ds[varname].sel(time=ds.time.dt.month.isin(7))[::3,:,:].data

def read_era5_grib(year,varname,unit,suffix):
    outfile = '%s/ERA5-%s_%d.nc'%(outdir,varname,year)
    if os.path.exists(outfile):
        print('%s exists'%outfile)
        return
    else:
        print('handle %s'%outfile)

    time = []
    var = []
    fn_stream = subprocess.check_output(
        'ls %s/%d0[78]*-%s.grib'%(
        indir[0],year,suffix), shell=True).decode('utf-8')
    fn_list = fn_stream.split()[0:36]
    
    for itm in fn_list:
        print(itm)
        ds = xr.open_dataset(itm,engine='cfgrib',backend_kwargs={
            'filter_by_keys': {'typeOfLevel': 'surface'}})
        time= time+ np.datetime_as_string(ds['time'].data,
            unit='h').tolist()
        var = var + ds[varname].data.tolist()
    var = np.array(var)
    time = pd.to_datetime(np.array(time),format='%Y-%m-%dT%H')
    print('var shape %s, time %s'%(str(var.shape), str(time.shape)))
    lat = ds['latitude'].data
    lon = ds['longitude'].data
    da = xr.DataArray(var.data,coords=[time,lat,lon],dims=['time','lat','lon'])
    da.attrs['units']=unit
    ds = da.to_dataset(name=varname)
    ds.to_netcdf(outfile,"w")

def calc_wrf_interp(year,varname,dh,path,casename):
    outfile = '%s/%s-%s_%d.nc'%(outdir,casename,varname,year)
    if os.path.exists(outfile):
        print('%s exists'%outfile)
        return
    else:
        print('handle %s'%outfile)

    fn_stream = subprocess.check_output(
        'ls %s%d070100/wrfout_d01_%d-*'%(
        path,year,year), shell=True).decode('utf-8')
    fn_list = fn_stream.split()
    
    wrf_list=[Dataset(itm) for itm in fn_list[::dh]]
    print('%d total file %d; opened file %d'%(year, 
        len(fn_list),len(wrf_list)))

    term = wrf.getvar(wrf_list, varname, timeidx=wrf.ALL_TIMES, 
        method='cat')
    time = wrf.getvar(wrf_list,'Times',timeidx=wrf.ALL_TIMES, 
        method='cat').data
    time = pd.to_datetime(np.datetime_as_string(
        time,unit='h'),format='%Y-%m-%dT%H')
    xlat = wrf.getvar(wrf_list[0],'XLAT').data
    xlon = wrf.getvar(wrf_list[0],'XLONG').data
    del wrf_list
    gc.collect()
    var = xr.apply_ufunc(grid_interp,term,xlat,xlon,
        input_core_dims=[['south_north','west_east'],
        ['south_north','west_east'],['south_north','west_east']],
        output_core_dims=[['lat','lon']],
        vectorize=True, dask="parallelized",output_dtypes=['float64'],)
                        
    da = xr.DataArray(var.data,coords=[time,ilat,ilon],dims=['time','lat','lon'])
    da.attrs['units']='K'
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
