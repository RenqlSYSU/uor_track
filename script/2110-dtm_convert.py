'''
read Whole_HK_DTM_5m.asc, 
Interpolate to a resolution of 100m
calculate latitude and longitude,
convert to nc file

20211011
'''

import xarray as xr
import pandas as pd
import numpy as np
import subprocess
from hk80 import LatLon, HK80

filname="/home/lzhenn/cooperate/data/Whole_HK_DTM_5m.asc"
fileout0="/home/lzhenn/cooperate/data/Whole_HK_DTM_5m.nc"
fileout1="/home/lzhenn/cooperate/data/Whole_HK_DTM_100m.nc"
data = np.loadtxt(filname, skiprows=6)
dim = data.shape # (nrows, ncols)
xll = 799997.5
yll = 799997.5
size0  = 5 #m
x0=np.arange(xll,xll+size0*dim[1],size0)
y0=np.arange(yll,yll+size0*dim[0],size0)
print(data.shape)
print("x numb: %d"%len(x0))
print("y numb: %d"%len(y0))

size1 = 100 #m
scale = int(size1/size0)
x1=np.arange(xll+size1/2.,x0[-1],size1)
y1=np.arange(yll+size1/2.,y0[-1],size1)
data1=np.empty((len(y1),len(x1)), dtype = float)
lat1=np.empty((len(y1)), dtype = float)
lon1=np.empty((len(x1)), dtype = float)
print("scale is %d"%scale)
print(data1.shape)
print("new x numb: %d"%len(x1))
print("new y numb: %d"%len(y1))

for nx in range(0,len(x1),1):
    hku = HK80(northing=y1[0], easting=x1[nx]).to_wgs84()
    lon1[nx] = hku.longitude
    for ny in range(0,len(y1),1):
        data1[ny,nx] = np.mean(data[scale*ny:scale*(ny+1),scale*nx:scale*(nx+1)])

for ny in range(0,len(y1),1):
    hku = HK80(northing=y1[ny], easting=x1[0]).to_wgs84()
    lat1[ny] = hku.latitude

ds = xr.Dataset(
         {
                 "dtm" : (["lat","lon"], data1),
                 },
         coords={
                 "lat" : (["lat"], lat1[::-1]),
                 "lon" : (["lon"], lon1),
                 },
         )
ds.attrs["description"] = "Space resolution %dm. "%size1+\
        "The data is computed from https://data.gov.hk/en-data/dataset/hk-landsd-openmap-5m-grid-dtm"
ds.to_netcdf(fileout1,"w")
del ds

lat0=np.empty((len(y0)), dtype = float)
lon0=np.empty((len(x0)), dtype = float)
for nx in range(0,len(x0),1):
    hku = HK80(northing=y0[0], easting=x0[nx]).to_wgs84()
    lon0[nx] = hku.longitude

for ny in range(0,len(y0),1):
    hku = HK80(northing=y0[ny], easting=x0[0]).to_wgs84()
    lat0[ny] = hku.latitude

ds = xr.Dataset(
         {
                 "dtm" : (["lat","lon"], data),
                 },
         coords={
                 "lat" : (["lat"], lat0[::-1]),
                 "lon" : (["lon"], lon0),
                 },
         )
ds.attrs["description"] = "Space resolution %dm. "%size0+\
        "The data is computed from https://data.gov.hk/en-data/dataset/hk-landsd-openmap-5m-grid-dtm"
ds.to_netcdf(fileout0,"w")

