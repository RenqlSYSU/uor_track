import matplotlib
matplotlib.use('Agg')

import numpy as np
import netCDF4
from datetime import datetime
import pyroms
import pyroms_toolbox
import sys

# get HYCOM grid data information 
invarname = 'water_temp'
outvarname = 'temp'

#read grid and variable attributes from the first file
url='/home/lzhenn/drv_field/hycom_subset/2021091500/hycom.grid.exp930.nc'
dataset = netCDF4.Dataset(url)

lon1 = dataset.variables['lon'][:]
lat1 = dataset.variables['lat'][:]
z = dataset.variables['depth'][:]
lon, lat = np.meshgrid(lon1, lat1)

var = dataset.variables[invarname][0,:,:,:]
spval = var.get_fill_value()
units = dataset.variables[invarname].units
long_name = dataset.variables[invarname].long_name

dataset.close()

year = 2007
day = 1

#create netCDF file
outfile = '/home/lzhenn/drv_field/hycom_subset/2021091500/hycom.njord.grid.exp930.nc'
nc = netCDF4.Dataset(outfile, 'w', format='NETCDF3_64BIT')
nc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
nc.title = 'HYCOM + NCODA Global 1/12 Analysis (GLBa0.08)'

#create dimensions
Mp, Lp = lon.shape 
N = len(z)
nc.createDimension('lon', Lp)
nc.createDimension('lat', Mp)
nc.createDimension('z', N)
nc.createDimension('ocean_time', None)

#create variables        
nc.createVariable('lon', 'f', ('lat', 'lon'))
nc.variables['lon'].long_name = 'longitude'
nc.variables['lon'].units = 'degrees_east'
nc.variables['lon'][:] = lon

nc.createVariable('lat', 'f', ('lat', 'lon'))
nc.variables['lat'].long_name = 'latitude'
nc.variables['lat'].units = 'degrees_north'
nc.variables['lat'][:] = lat

nc.createVariable('z', 'f', ('z'))
nc.variables['z'].long_name = 'depth'
nc.variables['z'].units = 'meter'
nc.variables['z'][:] = z

nc.createVariable('ocean_time', 'f', ('ocean_time'))
nc.variables['ocean_time'].units = 'days since 1900-01-01 00:00:00'
jday = pyroms_toolbox.date2jday(datetime(year, 1, 1)) + day - 1
nc.variables['ocean_time'][0] = jday

nc.createVariable(outvarname, 'f', ('ocean_time', 'z', 'lat', 'lon'), fill_value=spval)
nc.variables[outvarname].long_name = long_name
nc.variables[outvarname].units = units
nc.variables[outvarname][0] = var
        
nc.close()

