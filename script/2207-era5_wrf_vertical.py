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
years = range(2011,2021)
lats = 15
latn = 50
lonl = 110
lonr = 120
ilat = np.arange(lats, latn+0.1, 0.25)
ilon = np.arange(lonl, lonr+0.1, 0.25)
lev = [1000,925,850,700,600,500,400,300,250,200,100]

def main_run():
    #read_era5_grib('u','ua','m/s','pl')
    #for i in range(1,3):
    #    calc_wrf_interp('ua','m/s',3,indir[i],case[i])
    draw_clim_tcc('ua','U',np.arange(-20,31,3),np.arange(-2.4,2.5,0.3))

def draw_clim_tcc(varname,drawvar,cnlev1,cnlev2):
    var = [] 
    for ca in range(len(case)):
        ds = xr.open_dataset('%s/%s-%s_clim.nc'%(outdir,case[ca],varname))
        var.append(ds[varname].sel(lat=ilat).data)

    draw_3plot([var[1]-var[0],var[2]-var[0],var[2]-var[1]], 
        ['diff WRF-ERA5','diff WRFROMS-ERA5','diff WRFROMS-WRF'], 
        cnlev2, cmaps.BlueDarkRed18, '%s/cwrf_%s_diff.png'%(figdir,varname))
    draw_3plot([var[1]-var[0],var[2]-var[0],var[2]-var[1]], 
        ['diff WRF-ERA5','diff WRFROMS-ERA5','diff WRFROMS-WRF'], 
        np.arange(-0.8,0.81,0.1), cmaps.BlueDarkRed18, '%s/cwrf_%s_diff2.png'%(figdir,varname))
    draw_3plot(var, ['%s %s'%(drawvar,nc) for nc in case],cnlev1, 
        cmaps.precip2_17lev, '%s/cwrf_%s_clim.png'%(figdir,varname))

def draw_3plot(var,figtitle,cnlevels,fcolors,figfile):
    lat_sp = 10
    nrow = len(var) #2 #1 #6 #1 
    ncol = 2 
    norm = colors.BoundaryNorm(boundaries=cnlevels, 
        ncolors=fcolors.N,extend='both')
    
    fig = plt.figure(figsize=(12,12),dpi=300)
    ax = fig.subplots(nrow, ncol)
    for nr in range(0,len(var),1):
        print(var[nr].min())
        print(var[nr].max())
        axe = ax[nr][0]
        axe.set_title(figtitle[nr],fontsize=title_font)
        cont = axe.contourf(ilat, lev, var[nr], cnlevels, 
             cmap=fcolors,extend='both',norm=norm)
        plt.colorbar(cont, ax=axe, shrink=.9, pad=0.01)
        
        axe.set_ylabel('lev(hPa)',fontsize=label_font,fontdict=font)
        axe.set_yscale('symlog')
        axe.set_yticklabels(np.arange(1000,5,-100))
        axe.set_yticks(np.arange(1000,5,-100))
        axe.set_ylim(1000,100)
        
        axe.set_xlabel('',fontsize=label_font,fontdict=font)
        axe.set_xticks(np.arange(lats,latn,lat_sp))
        axe.xaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
    
    plt.savefig(figfile, bbox_inches='tight',pad_inches=0.01)
    
def read_era5_grib(varname,new_var,unit,suffix):
    outfile = '%s/ERA5-%s_clim.nc'%(outdir,new_var)
    if os.path.exists(outfile):
        print('%s exists'%outfile)
        return
    else:
        print('handle %s'%outfile)
    
    numb = 0
    for year in years:
        fn_stream = subprocess.check_output(
            'ls %s/%d07*-%s.grib'%(
            indir[0],year,suffix), shell=True).decode('utf-8')
        fn_list = fn_stream.split()
        for itm in range(len(fn_list)):
            print(fn_list[itm])
            ds = xr.open_dataset(fn_list[itm],engine='cfgrib')
            if year==2011 and itm==0:
                var = ds[varname].sel(isobaricInhPa=lev,
                    longitude=ilon).mean(('time','longitude'))
            else:
                var.data = var.data + ds[varname].sel(isobaricInhPa=lev,
                    longitude=ilon).mean(('time','longitude')).data
            numb = numb + 1
    print('var shape %s'%(str(var.shape)))
    da = xr.DataArray(var.data/numb,coords=[lev,var.latitude],
        dims=['lev','lat'])
    da.attrs['units']=unit
    ds = da.to_dataset(name=new_var)
    ds.to_netcdf(outfile,"w")

def calc_wrf_interp(varname,unit,dh,path,casename):
    outfile = '%s/%s-%s_clim.nc'%(outdir,casename,varname)
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
        p = wrf.getvar(wrf_list, 'pressure', timeidx=wrf.ALL_TIMES, 
            method='cat')
        if year == 2011:
            term = wrf.interplevel(z, p, lev).mean('Time')
            xlat = wrf.getvar(wrf_list[0],'XLAT')
            xlon = wrf.getvar(wrf_list[0],'XLONG')
            print(term)
        else:
            term.data = term.data + wrf.interplevel(z, p, lev).mean('Time').data
        del wrf_list
        gc.collect()

    term.data = term.data/len(years)
    var = xr.apply_ufunc(grid_interp,term,xlat,xlon,
        input_core_dims=[['south_north','west_east'],
        ['south_north','west_east'],['south_north','west_east']],
        output_core_dims=[['lat','lon']],
        vectorize=True, dask="parallelized",output_dtypes=['float64'],)
    print(var)
                        
    da = xr.DataArray(var.data.mean(axis=2),
        coords=[lev,ilat],dims=['lev','lat'])
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
