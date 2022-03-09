import xarray as xr
import pandas as pd
import numpy as np
from scipy import interpolate
from datetime import datetime
import subprocess

def inter2d(var, dst_lat, dst_lon):
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

def inter_depth(var,z_rho):
    dst = var.interp(depth=(-1*z_rho),method='linear')
    return dst.data

def average_depth(var,z_rho):
    dst = (var*np.diff(z_w)).sum()/(-1*z_w[0])
    return dst

def inter3d(var, dst_lat, dst_lon):
    zeta = ds_roms.zeta[0,:,:]
    h = zeta
    h.data = ds_grid.h.data
    z_rho = sigma2depth(zeta,h,ds_roms.sc_r,ds_roms.Cs_r)
    
    var = var.interpolate_na(dim="depth", 
        method="linear", fill_value="extrapolate")
    var = var.interpolate_na(dim="lon", 
        method="linear", fill_value="extrapolate")
    
    var = var.interp(lat=dst_lat,lon=dst_lon,
        method='linear').data 
    #dst_value = xr.apply_ufunc(inter_depth,var,z_rho,
    #    input_core_dims=[['depth'], ['sc_r']],
    #    vectorize=True,
    #    dask="parallelized",
    #    output_dtypes=['float64'],)
    dst_value = ds_roms['temp'].data
    for i in range(len(dst_lat)):
        for j in range(len(dst_lon)):
            dst_value[0,:,i,j] = var[0,:,i,j].interp(
                depth=(-1*z_rho[i,j,:]),method='linear')
    return dst_value

def inter_uv(var,dst_lat,dst_lon,varname):
    term = ds_hycom['surf_el'][0,:,:].interpolate_na(
        dim="lon", method="linear",fill_value="extrapolate")
    zeta = ds_roms[varname][0,0,:,:]
    zeta.data = term.interp(lat=dst_lat,lon=dst_lon,
        method='linear').data 
    x_rho = ds_grid['x_rho'].data
    y_rho = ds_grid['y_rho'].data
    h_rho = ds_grid['h_rho'].data
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
        method='linear').data 
    dst_value = ds_roms[varname].data
    dst_bar = ds_roms[varname+'bar'].data
    for i in range(len(dst_lat)):
        for j in range(len(dst_lon)):
            dst_value[0,:,i,j] = var[0,:,i,j].interp(
                depth=(-1*z_rho[i,j,:]),method='linear')
            dst_bar[0,i,j] = (dst_value[0,:,i,j] * 
                np.diff(z_w[i,j,:])).sum() / (-1*z_w[i,j,0])
    return dst_value, dst_bar

dst_dir='/home/lzhenn/cooperate/script/roms_drv/'
out_file  = 'coawst_clm_20210915.nc' 

hycom_dir = '/home/lzhenn/drv_field/hycom_subset/2021091500/'
roms_dir  = '/home/lzhenn/drv_field/icbc/2021091500/'
roms_grid = '/home/lzhenn/Njord/Projects/Njord/roms_swan_grid/roms_d01_lp0d1.nc'

ds_hycom = xr.open_dataset('%shycom_glby_930_2021091500.nc'%hycom_dir)
ds_roms  = xr.open_dataset('%scoawst_ini.nc'%roms_dir)
ds_grid  = xr.open_dataset(roms_grid)
dst_ds   = xr.open_dataset('%s%s'%(roms_dir,out_file))

#ds_roms['ocean_time'].data = ds_hycom

lat_rho = ds_grid['lat_rho'][:,0].data
lon_rho = ds_grid['lon_rho'][0,:].data
print(ds_roms['zeta'][0,30:32,40:42])
term = ds_hycom['surf_el'].interpolate_na(
    dim="lon", method="linear",fill_value="extrapolate")
term1 = term.interp(lat=lat_rho,lon=lon_rho,method='linear').data 
term2 = ds_roms['zeta'].data
print((term1-term2)/term2)
ds_roms['zeta'].data = term1
#ds_roms['zeta'].data = term.interp(
#   lat=lat_rho,lon=lon_rho,method='linear').data 
print(ds_roms['zeta'][0,30:32,40:42])

term = ds_hycom['water_temp'][0,:,:,:].interpolate_na(
        dim="depth", method="linear",fill_value="extrapolate")
term = term.interpolate_na(
    dim="lon", method="linear",fill_value="extrapolate")
term = term.interp(lat=lat_rho,lon=lon_rho,method='linear')
print(term)

zeta = ds_roms.zeta[0,:,:]
h = zeta
h.data = ds_grid.h.data
if ds_roms.Vtransform == 1:
    Zo_rho = ds_roms.hc.data * (ds_roms.sc_r - ds_roms.Cs_r) + ds_roms.Cs_r * h 
    z_rho = Zo_rho + zeta * (1 + Zo_rho / h)
elif ds_roms.Vtransform == 2:
    Zo_rho = (ds_roms.hc.data * ds_roms.sc_r + ds_roms.Cs_r * h) / (ds_roms.hc.data + h)
    z_rho = zeta + Zo_rho * (zeta + h)

print(ds_roms['temp'][0,:,30:32,40:42].data)
for i in range(len(lat_rho)):
    for j in range(len(lon_rho)):
        var = term[:,i,j].interp(
            depth=(-1*z_rho[i,j,:]),method='linear')
        ds_roms['temp'][0,:,i,j].data = var
print(ds_roms['temp'][0,:,30:32,40:42].data)

ds_roms.to_netcdf('%s%s'%(dst_dir,out_file),"w")

