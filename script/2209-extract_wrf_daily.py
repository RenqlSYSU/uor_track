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

paths = ['/home/metctm1/array_hq133/data/s2s/wrfonly/realtime',
         '/home/metctm1/array_hq133/data/s2s/wrfroms/realtime']
indir = ['%s/WRF_NG'%paths[0],
         '%s/WRF_free'%paths[0],
         '%s/CFS_WRF_free'%paths[0],
         '%s/CP_NG'%paths[1],
         '%s/CP_free'%paths[1],
         '%s/CFS_CP_free'%paths[1],
         '%s/CFSA_HYCOM'%paths[1],
         '%s/ERA5_CFSO'%paths[1]]
case = [x.split('/')[-1] for x in indir] 
years = [2022,]#
#indir = ['/home/metctm1/array_hq133/data/s2s/wrfonly/clim/',
#         '/home/metctm1/array_hq133/data/s2s/wrfroms/clim/aegir_']
#case = ['WRF','WRFROMS']
#years = range(2011,2021) #[2022,]#
outdir = '/home/lzhenn/cooperate/data/cwrf2'
figdir = '/home/lzhenn/cooperate/fig'
ftime  = pd.date_range(start='2021-07-01 00',end='2021-08-04 23',
    freq='1D',closed=None)
print('%d days'%len(ftime))

def main_run():
    for year in years:
        for i in range(len(case)):
            calc_wrf_interp('T2','K',3,indir[i],case[i],year)
            calc_wrf_precip('tp','mm/day',indir[i],case[i],year)

def calc_wrf_precip(varname,unit,path,casename,year):
    outfile = '%s/%s-%s_daily_%d.nc'%(outdir,casename,varname,year)
    if os.path.exists(outfile):
        print('%s exists'%outfile)
        return
    else:
        print('handle %s'%outfile)
    
    fn_stream = subprocess.check_output(
        'ls %s/wrfout_d01_%d-0[78]-*_00:00:00'%(
        path,year), shell=True).decode('utf-8')
    #fn_stream = subprocess.check_output(
    #    'ls %s%d070100/wrfout_d01_%d-*_00:00:00'%(
    #    path,year,year), shell=True).decode('utf-8')
    fn_list = fn_stream.split()
    wrf_list=[Dataset(itm) for itm in fn_list]
    print('%d total file %d; opened file %d'%(year, 
        len(fn_list),len(wrf_list)))
    
    z = wrf.getvar(wrf_list, 'RAINNC', timeidx=wrf.ALL_TIMES, 
        method='cat')
    z.data = wrf.getvar(wrf_list, 'RAINC', timeidx=wrf.ALL_TIMES, 
        method='cat').data + z.data
    term = z[:-1]
    term.data = z.data[1:]-z.data[:-1]
    term.attrs['projection']='LambertConformal'
    term.attrs['units']=unit
    print(term)
    del wrf_list, z
    gc.collect()

    ds = term.to_dataset(name=varname)
    ds.to_netcdf(outfile,"w")

def calc_wrf_interp(varname,unit,dh,path,casename,year):
    outfile = '%s/%s-%s_daily_%d.nc'%(outdir,casename,varname,year)
    if os.path.exists(outfile):
        print('%s exists'%outfile)
        return
    else:
        print('handle %s'%outfile)
    
    fn_stream = subprocess.check_output(
        'ls %s/wrfout_d01_%d-0[78]-*'%(
        path,year), shell=True).decode('utf-8')
    #fn_stream = subprocess.check_output(
    #    'ls %s%d070100/wrfout_d01_%d-*'%(
    #    path,year,year), shell=True).decode('utf-8')
    fn_list = fn_stream.split()
    wrf_list=[Dataset(itm) for itm in fn_list[:-1:dh]]
    print('%d total file %d; opened file %d'%(year, 
        len(fn_list),len(wrf_list)))
    
    z = wrf.getvar(wrf_list, varname, timeidx=wrf.ALL_TIMES, 
        method='cat')
    term = z.groupby(z.Time.dt.dayofyear).mean('Time').rename({'dayofyear':'time'})
    term.coords['time'] = ftime
    print(term)
    term.attrs['projection']='LambertConformal'
    term.attrs['units']=unit
    del wrf_list, z
    gc.collect()

    ds = term.to_dataset(name=varname)
    ds.to_netcdf(outfile,"w")

if __name__=='__main__':
    main_run()

