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
    #for year in range(2011,2021):
    #    read_era5_grib(year,'t2m','T2','K','sl')
    #    read_era5_grib(year,'skt','TSK','K','sl')
    #    read_era5_grib(year,'z','500z','m2*s-2','pl')
        #read_era5_grib(year,'u','200u','m/s','pl')
    #    for i in range(1,3):
            #calc_wrf_interp(year,'T2',3,indir[i],case[i])
    #        calc_wrf_interp(year,'TSK',3,indir[i],case[i])
            #calc_wrf_interp(year,'500z',3,indir[i],case[i])
     #       calc_wrf_interp(year,'200u',3,indir[i],case[i])
     #       calc_wrf_interp_rainfall(year,'tp',744,indir[i],case[i])
    #draw_clim_tcc('TSK','TSK',np.arange(0,34,2),np.arange(-4,4.1,0.5))
    #draw_clim_tcc('T2','T2m',np.arange(0,34,2),np.arange(-4,4.1,0.5))
    draw_clim_tcc('500z','500hPa z',np.arange(5740,5910,10),np.arange(-8,8.1,1))
    #draw_clim_tcc('200u','200hPa U',np.arange(-4,30,2),np.arange(-4,4.1,0.5))
    #draw_clim_tcc('tp','Tp(mm/day)',np.arange(0,34,2),np.arange(-8,8.1,1))

def draw_clim_tcc(varname,drawvar,cnlev1,cnlev2):
    var = [] 
    for ca in range(len(case)):
        ds = xr.open_mfdataset('%s/%s-%s_20*.nc'%(outdir,case[ca],varname),
           combine='nested', concat_dim='time', parallel=True)
        if varname == 'tp' and ca == 0:
            var.append(ds[varname].sel(latitude=ilat,longitude=ilon).load())
            #var[0] = var[0].drop('time',dim=None)
            var[0].data = var[0].data*31
        else:
            var.append(ds[varname].sel(lat=ilat,lon=ilon).load())
        #var[ca].data = np.nan_to_num(var[ca].data, nan=0)
        #var.append(ds[varname].sel(lat=ilat,lon=ilon).data.mean(axis=0)-273.15)
    ''' 
    region = regionmask.defined_regions.natural_earth_v5_0_0.countries_110
    mask = region.mask(var[0]).data
    var = [x.where(mask==139) for x in var]
    '''
    #shf = '/disk/r074/lzhenn/tracacode/2205-PPOL-COOP/data/gadm36_CHN_0.shp'
    shf = '/home/metctm1/array/tracacode/UTILITY-2016/shp/cnmap/gadm36_CHN_1.shp'
    mask = regionmask.mask_geopandas(gpd.read_file(shf),var[0])
    var = [x.where(mask==5) for x in var]
    print(mask)
    print(var[0])

    if varname == 'tp':
        term1 = var[0].data
        var[0] = var[1].copy()
        var[0].data = term1
    if varname in ['500z','200u']:
        var[0] = var[0]/9.8#.data
        var[1] = var[1]*9.8#.data
        var[2] = var[2]*9.8#.data
    '''
    term = []
    term.append(xr.apply_ufunc(calc_corr,var[0],var[1],
        input_core_dims=[['time'],['time']],
        vectorize=True, dask="parallelized",output_dtypes=['float64'],))
    term.append(xr.apply_ufunc(calc_corr,var[0],var[2],
        input_core_dims=[['time'],['time']],
        vectorize=True, dask="parallelized",output_dtypes=['float64'],))
    term.append(xr.apply_ufunc(calc_corr,var[2],var[1],
        input_core_dims=[['time'],['time']],
        vectorize=True, dask="parallelized",output_dtypes=['float64'],))
    draw_3plot(term, ['TCC WRF & ERA5','TCC WRFROMS & ERA5','TCC WRFROMS & WRF'], 
        np.concatenate((np.arange(-0.95,-0.56,0.05),np.array([0]),np.arange(0.6,0.96,0.05))), 
        cmaps.BlueDarkRed18, '%s/cwrf_%s_corr.png'%(figdir,varname))
    del term
    '''
    if varname == 'tp':
        var = [x.data.mean(axis=0)*1000/31.0 for x in var]
    elif varname in ['T2','TSK']:
        var = [x.data.mean(axis=0)-273.15 for x in var]
    else:
        var = [x.mean(axis=0) for x in var]

#    draw_3plot([var[1]-var[0],var[2]-var[0],var[2]-var[1]], 
#        ['diff WRF-ERA5','diff WRFROMS-ERA5','diff WRFROMS-WRF'], 
#        cnlev2, cmaps.BlueDarkRed18, '%s/cwrf_%s_diff.png'%(figdir,varname))
#    draw_3plot([var[1]-var[0],var[2]-var[0],var[2]-var[1]], 
#        ['diff WRF-ERA5','diff WRFROMS-ERA5','diff WRFROMS-WRF'], 
#        np.arange(-0.8,0.81,0.1), cmaps.BlueDarkRed18, '%s/cwrf_%s_diff2.png'%(figdir,varname))
#    
    draw_3plot(var, ['%s %s'%(drawvar,nc) for nc in case],cnlev1, 
        cmaps.precip2_17lev, '%s/cwrf_%s_clim.png'%(figdir,varname))

def draw_3plot(var,figtitle,cnlevels,fcolors,figfile):
    lat_sp = 15
    lon_sp = 20 #60 #
    nrow = len(var) #2 #
    ncol = 1 #6 #
    norm = colors.BoundaryNorm(boundaries=cnlevels, 
        ncolors=fcolors.N,extend='both')
    
    fig = plt.figure(figsize=(12,12),dpi=300)
    ax = fig.subplots(nrow, ncol, subplot_kw=dict(
        projection=ccrs.PlateCarree(central_longitude=180.0)))
    for nr in range(0,len(var),1):
        print(var[nr].min())
        print(var[nr].max())
        axe = ax[nr]
        axe.set_title(figtitle[nr],fontsize=title_font)
        
        axe.add_feature(cfeat.NaturalEarthFeature('physical', 'land', '50m',
            edgecolor='k',facecolor='none',linewidth=0.8))
        coast_shp2 = Reader(os.getenv("SHP_LIB")+"/cnmap/cnhimap.dbf").geometries()
        coastline2 = cfeat.ShapelyFeature(coast_shp2, ccrs.PlateCarree(), 
            edgecolor="k",facecolor="none")
        axe.add_feature(coastline2, linewidth=0.8,zorder=2)

        cont = axe.contourf(ilon, ilat, var[nr], cnlevels, 
             transform=ccrs.PlateCarree(),cmap=fcolors,
             extend='both',norm=norm)
        plt.colorbar(cont, ax=axe, shrink=.9, pad=0.01)
        
        axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
        axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
        axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
        axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))

    plt.savefig(figfile, bbox_inches='tight',pad_inches=0.01)
    
def read_era5_grib(year,varname,new_var,unit,suffix):
    outfile = '%s/ERA5-%s_%d.nc'%(outdir,new_var,year)
    if os.path.exists(outfile):
        print('%s exists'%outfile)
        return
    else:
        print('handle %s'%outfile)

    var = []
    fn_stream = subprocess.check_output(
        'ls %s/%d07*-%s.grib'%(
        indir[0],year,suffix), shell=True).decode('utf-8')
    fn_list = fn_stream.split()
    
    for itm in fn_list:
        print(itm)
        if new_var in ['T2','TSK']:
            ds = xr.open_dataset(itm,engine='cfgrib',backend_kwargs={
                'filter_by_keys': {'typeOfLevel': 'surface'}})
            var = var + ds[varname].data.tolist()
        if new_var == '200u':
            ds = xr.open_dataset(itm,engine='cfgrib')
            var = var + ds[varname].sel(isobaricInhPa=200).data.tolist()
        if new_var == '500z':
            ds = xr.open_dataset(itm,engine='cfgrib')
            var = var + ds[varname].sel(isobaricInhPa=500).data.tolist()
    var = np.array(var).mean(axis=0)
    print('var shape %s'%(str(var.shape)))
    lat = ds['latitude'].data
    lon = ds['longitude'].data
    da = xr.DataArray(var.data,coords=[lat,lon],dims=['lat','lon'])
    da.attrs['units']=unit
    ds = da.to_dataset(name=new_var)
    ds.to_netcdf(outfile,"w")

def calc_wrf_interp(year,varname,dh,path,casename):
    outfile = '%s/%s-%s_%d.nc'%(outdir,casename,varname,year)
    if os.path.exists(outfile):
        print('%s exists'%outfile)
        return
    else:
        print('handle %s'%outfile)
    
    fn_stream = subprocess.check_output(
        'ls %s%d070100/wrfout_d01_%d-07-*'%(
        path,year,year), shell=True).decode('utf-8')
    fn_list = fn_stream.split()
    
    wrf_list=[Dataset(itm) for itm in fn_list[::dh]]
    print('%d total file %d; opened file %d'%(year, 
        len(fn_list),len(wrf_list)))
    
    if varname in ['T2','TSK']:
        term = wrf.getvar(wrf_list, varname, timeidx=wrf.ALL_TIMES, 
            method='cat').mean('Time').data
        unit = 'k'
    if varname == '200u':
        z = wrf.getvar(wrf_list, 'ua', timeidx=wrf.ALL_TIMES, 
            method='cat')
        p = wrf.getvar(wrf_list, 'pressure', timeidx=wrf.ALL_TIMES, 
            method='cat')
        term = wrf.interplevel(z, p, 200).mean('Time').data
        unit = 'm/s'
    if varname == '500z':
        z = wrf.getvar(wrf_list, 'z', timeidx=wrf.ALL_TIMES, 
            method='cat')
        p = wrf.getvar(wrf_list, 'pressure', timeidx=wrf.ALL_TIMES, 
            method='cat')
        term = wrf.interplevel(z, p, 500).mean('Time').data/9.8
        unit = 'gpm'
    xlat = wrf.getvar(wrf_list[0],'XLAT').data
    xlon = wrf.getvar(wrf_list[0],'XLONG').data
    del wrf_list
    gc.collect()
    var = grid_interp(term,xlat,xlon)
                        
    da = xr.DataArray(var.data,coords=[ilat,ilon],dims=['lat','lon'])
    da.attrs['units'] = unit
    ds = da.to_dataset(name=varname)
    ds.to_netcdf(outfile,"w")

def calc_wrf_interp_rainfall(year,varname,dh,path,casename):
    outfile = '%s/%s-%s_%d.nc'%(outdir,casename,varname,year)
    if os.path.exists(outfile):
        print('%s exists'%outfile)
        return
    else:
        print('handle %s'%outfile)
    itm = '%s%d070100/wrfout_d01_%d-08-01_00:00:00'%(path,year,year)
    wrf_list = Dataset(itm) 
    term = wrf.getvar(wrf_list,'RAINC').data+wrf.getvar(wrf_list,'RAINNC').data
    xlat = wrf.getvar(wrf_list,'XLAT').data
    xlon = wrf.getvar(wrf_list,'XLONG').data
    del wrf_list
    gc.collect()
    var = grid_interp(term,xlat,xlon)/1000
    da = xr.DataArray(var.data,coords=[ilat,ilon],dims=['lat','lon'])
    da.attrs['units'] = 'm/month'
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

def calc_corr(x,y):
    nas = np.logical_or(np.isnan(x), np.isnan(y))
    r,p = stats.pearsonr(x[~nas], y[~nas])
    return r

if __name__=='__main__':
    main_run()
