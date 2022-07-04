#!/usr/bin/env python
import xarray as xr
import numpy as np

lonl=70  #0  #
lonr=150#360#
lats=15  #
latn=60 #

fname = '/gws/nopw/j04/ncas_generic/users/renql/ERA5_mon/ERA5_mon_sp_1979-2020.nc'
ds  = xr.open_dataset(fname)
lat = ds.latitude
lon = ds.longitude
ilon = lon[(lon>=lonl) & (lon<=lonr)]
ilat = lat[(lat>=lats) & (lat<=latn)]
var = ds['sp'].sel(latitude=ilat,longitude=ilon).mean('time')

ds2 = var.to_dataset(name='sp')
ds2.to_netcdf("/home/users/qd201969/uor_track/mdata/ERA5_annual_sp_1979-2020.nc","w")
