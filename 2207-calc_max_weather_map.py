#!/usr/bin/env python
import xarray as xr
import numpy as np
import pandas as pd
from multiprocessing import Pool
import sys, os, subprocess, linecache, gc
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmaps
from renql import composite

prefix = "fftadd"
radiu = 2.5
title= {'_local':'local',
        '_outside':'outside',
        '_total':'total',
        '':'All'}
suffix = '_total'
lev  = [850,500,250]
path = '/home/users/qd201969/ERA5-1HR-lev'
datapath ="/home/users/qd201969/uor_track/mdata"
figdir = '/home/users/qd201969/uor_track/fig'

lonl=0  #0  #
lonr=150#360#
lats=15  #
latn=70 #

def main_run():
    #gene_grid()
    calc_monthly_max_mean(850,9,'%s/max_mean_10mwind_%d.nc'%(datapath,850))

def calc_monthly_max_mean(nl,nint,fileout):
    # nint, intensity columns in fftadd file
    # nint=9,max10mwind; nint=10,average precip
    ds = xr.open_dataset("%s/grid_%.1f.nc"%(datapath,radiu))
    grid = ds['grid']
    dim = grid.shape
    print(dim)
    mean = np.zeros([12,dim[0],dim[1]],dtype=float)
    maxv = np.zeros([12,dim[0],dim[1]],dtype=float)
    numb = np.zeros([12,dim[0],dim[1]],dtype=float)

    ctime,clat,clon,cinte = composite.composite_time(
        "%s/%s_%d_1980-2020%s"%(
        path,prefix,nl,suffix),lats,latn,lonl,lonr,nint)
    nd = 0
    for ct,clo,cla,cit in zip(ctime,clon,clat,cinte):
        nd = nd + 1
        dist = np.square(grid.lat-cla)+np.square(grid.lon-clo)
        ind = np.unravel_index(np.argmin(dist.data), dist.shape)
        print('%d dist shape %s, ind %s'%(nd,str(dist.shape),str(ind)))
        numb[ct.month-1,ind[0],ind[1]] = numb[ct.month-1,ind[0],ind[1]] + 1
        mean[ct.month-1,ind[0],ind[1]] = mean[ct.month-1,ind[0],ind[1]] + cit
        if cit > maxv[ct.month-1,ind[0],ind[1]]:
            maxv[ct.month-1,ind[0],ind[1]] = cit

    mean = mean/numb
    ds = xr.Dataset(
            {
                "numb": (["month", "lat", "lon"], numb),
                "maxv": (["month", "lat", "lon"], maxv),
                "mean": (["month", "lat", "lon"], mean),
                },
            coords={
                "month": range(0,12,1), 
                "lat"  : (["lat"],grid.lat),
                "lon"  : (["lon"],grid.lon),
                },
            )
    ds.attrs["description"]='%d cyclone point'%(len(ctime))
    ds.to_netcdf(fileout,"w")

def gene_grid():
    lat = np.arange(lats,latn+0.1,radiu)
    lon = np.arange(lonl,lonr+0.1,radiu)
    var = np.zeros([len(lat),len(lon)],dtype=float)
    da = xr.DataArray(var, coords=[lat,lon], dims=["lat", "lon"])
    ds = da.to_dataset(name="grid")
    ds.to_netcdf("%s/grid_%.1f.nc"%(datapath,radiu),"w")

if __name__=='__main__':
    main_run()
