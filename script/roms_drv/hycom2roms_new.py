import xarray as xr
import pandas as pd
import numpy as np
from scipy import interpolate
from datetime import datetime
import subprocess

def inter2d(var, dst_lat, dst_lon):
    '''
    first fill the missing value
    then interpolate to new 2d grid 
    '''
    print("interpolate %s"%var.long_name)
    var = var.interpolate_na(dim="lon", 
        method="linear", fill_value="extrapolate")
    dst_value = var.interp(lat=dst_lat,lon=dst_lon,
        method='linear').data 
    return dst_value

def sigma2depth(zeta,h,sc,Cs):
    '''
    convert S-coordinate to depth (m) 
    zeta: 2D free-surface
    h: 2D Final bathymetry, same dimension as zeta
    sc: 1D S-coordinate at RHO-points or W-points
    Cs: 1D S-coordinate stretching curves 
        at RHO-points or W-points
    return z_rho(erho: 514, xrho: 770, sc_r: 30)
    '''
    hc = ds_roms.hc.data
    # S-coordinate parameter, critical depth
    if ds_roms.Vtransform == 1:
        Zo_rho = hc * (sc - Cs) + Cs * h 
        z_rho = Zo_rho + zeta * (1 + Zo_rho / h)
    elif ds_roms.Vtransform == 2:
        Zo_rho = (hc * sc + Cs * h) / (hc + h)
        z_rho = zeta + Zo_rho * (zeta + h)
    return z_rho

def inter3d(var, dst_lat, dst_lon, ds_hycom):
    print("interplate %s"%var.long_name)
    term = ds_hycom['surf_el'][0,:,:].interpolate_na(
        dim="lon", method="linear",fill_value="extrapolate")
    zeta = ds_roms[varname][0,0,:,:]
    zeta.data = term.interp(lat=dst_lat,lon=dst_lon,
        method='linear').data 
    h = zeta
    h.data = ds_grid.h.data
    z_rho = sigma2depth(zeta,h,ds_roms.sc_r,ds_roms.Cs_r)
    
    var = var.interpolate_na(dim="depth", 
        method="linear", fill_value="extrapolate")
    var = var.interpolate_na(dim="lon", 
        method="linear", fill_value="extrapolate")
    
    var = var.interp(lat=dst_lat,lon=dst_lon,
        method='linear')
    dst_value = ds_roms['temp'].data
    for i in range(len(dst_lat)):
        for j in range(len(dst_lon)):
            dst_value[0,:,i,j] = var[0,:,i,j].interp(
                depth=(-1*z_rho[i,j,:]),method='linear')
    return dst_value

def inter_uv(var, dst_lat, dst_lon, varname, ds_hycom):
    print("interplate %s"%var.long_name)
    term = ds_hycom['surf_el'][0,:,:].interpolate_na(
        dim="lon", method="linear",fill_value="extrapolate")
    zeta = ds_roms[varname][0,0,:,:]
    zeta.data = term.interp(lat=dst_lat,lon=dst_lon,
        method='linear').data 
    
    x_rho = ds_grid['x_rho'].data
    y_rho = ds_grid['y_rho'].data
    h_rho = ds_grid['h'].data
    f = interpolate.interp2d(x_rho, y_rho, h_rho, kind='cubic')
    x_u = ds_grid['x_'+varname].data
    y_u = ds_grid['y_'+varname].data
    h = zeta
    h.data = f(x_u,y_u)
    
    z_rho = sigma2depth(zeta,h,ds_roms.sc_r,ds_roms.Cs_r)
    z_w   = sigma2depth(zeta,h,ds_roms.sc_w,ds_roms.Cs_w)
    
    var = var.interpolate_na(dim="depth", 
        method="linear", fill_value="extrapolate")
    var = var.interpolate_na(dim="lon", 
        method="linear", fill_value="extrapolate")
    
    var = var.interp(lat=dst_lat,lon=dst_lon,
        method='linear')
    dst_value = ds_roms[varname].data
    dst_bar = ds_roms[varname+'bar'].data
    for i in range(len(dst_lat)):
        for j in range(len(dst_lon)):
            dst_value[0,:,i,j] = var[0,:,i,j].interp(
                depth=(-1*z_rho[i,j,:]),method='linear')
            dst_bar[0,i,j] = (dst_value[0,:,i,j] * 
                np.diff(z_w[i,j,:])).sum() / (-1*z_w[i,j,0])
    return dst_value, dst_bar

def make_ini(time):
    ds_hycom = xr.open_dataset('%shycom_glby_930_%s00.nc'%(hycom_dir,time))
    dst_ds = xr.open_dataset('%scoawst_ini.nc.nc'%roms_dir)
    
    #dst_ds['ocean_time'].data = ds_hycom['time']
    
    lat_rho = ds_grid['lat_rho'][:,0].data
    lon_rho = ds_grid['lon_rho'][0,:].data
    dst_ds['zeta'].data = inter2d(
        ds_hycom['surf_el'], lat_rho, lon_rho)
    dst_ds['temp'].data = inter3d(
        ds_hycom['water_temp'], lat_rho, lon_rho, ds_hycom)
    dst_ds['salt'].data = inter3d(
        ds_hycom['salinity'], lat_rho, lon_rho, ds_hycom)

    lat_rho = ds_grid['lat_u'][:,0].data
    lon_rho = ds_grid['lon_u'][0,:].data
    dst_ds['u'].data, dst_ds['ubar'].data = inter_uv(
        ds_hycom['water_u'], lat_rho, lon_rho, 'u', ds_hycom)

    lat_rho = ds_grid['lat_v'][:,0].data
    lon_rho = ds_grid['lon_v'][0,:].data
    dst_ds['v'].data, dst_ds['vbar'].data = inter_uv(
        ds_hycom['water_v'], lat_rho, lon_rho, 'v', ds_hycom)

    dst_ds.to_netcdf('%scoawst_ini.nc'%dst_dir,"w")

def make_clm(time):
    ds_hycom = xr.open_dataset('%shycom_glby_930_%s00.nc'%(hycom_dir,time))
    dst_ds = xr.open_dataset('%scoawst_clm_20210915.nc'%roms_dir)
    
    dst_ds.to_netcdf('%scoawst_clm_%s.nc'%(dst_dir,time),"w")

def make_bdry(time):
    ds_hycom = xr.open_dataset('%shycom_glby_930_%s00.nc'%(hycom_dir,time))
    dst_ds = xr.open_dataset('%scoawst_bdy_20210915.nc'%roms_dir)
    
    dst_ds.to_netcdf('%scoawst_bdy_%s.nc'%(dst_dir,time),"w")

time = '20210915'
dst_dir='/home/lzhenn/cooperate/script/roms_drv/'
hycom_dir = '/home/lzhenn/drv_field/hycom_subset/2021091500/'

# fixed reference roms file
roms_dir  = '/home/lzhenn/drv_field/icbc/2021091500/'
roms_grid = '/home/lzhenn/Njord/Projects/Njord/roms_swan_grid/roms_d01_lp0d1.nc'
ds_roms  = xr.open_dataset('%scoawst_ini.nc'%roms_dir)
ds_grid  = xr.open_dataset(roms_grid)

lat_rho = ds_grid['lat_rho'][:,0].data
lon_rho = ds_grid['lon_rho'][0,:].data
term1 = ds_roms['zeta'].data
ds_roms['zeta'].data = inter2d(
    ds_hycom['surf_el'], lat_rho, lon_rho)
term2 = ds_roms['zeta'].data
print('difference percent:')
print((term2-term1)/term1)

term1 = ds_roms['temp'].data
ds_roms['temp'].data = inter3d(
    ds_hycom['water_temp'], lat_rho, lon_rho)
term2 = ds_roms['temp'].data
print('difference percent:')
print((term2-term1)/term1)

ds_roms['salt'].data = inter3d(
    ds_hycom['salinity'], lat_rho, lon_rho)

term1 = ds_roms['u'].data
lat_rho = ds_grid['lat_u'][:,0].data
lon_rho = ds_grid['lon_u'][0,:].data
ds_roms['u'].data, ds_roms['ubar'].data = inter_uv(
    ds_hycom['water_u'], lat_rho, lon_rho, 'u')
term2 = ds_roms['u'].data
print('difference percent:')
print((term2-term1)/term1)

term1 = ds_roms['vbar'].data
lat_rho = ds_grid['lat_v'][:,0].data
lon_rho = ds_grid['lon_v'][0,:].data
ds_roms['v'].data, ds_roms['vbar'].data = inter_uv(
    ds_hycom['water_v'], lat_rho, lon_rho, 'v')
term2 = ds_roms['vbar'].data
print('difference percent:')
print((term2-term1)/term1)

ds_roms.to_netcdf('%s%s'%(dst_dir,out_file),"w")

