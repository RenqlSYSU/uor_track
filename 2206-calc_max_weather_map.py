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

prefix = "fftadd"
radiu = 2.5
title= {'_local':'local',
        '_outside':'outside',
        '_total':'total',
        '':'All'}

fileout="/home/users/qd201969/uor_track/mdata"
lev  = [850,500,250]
path = '/home/users/qd201969/ERA5-1HR-lev'
figdir = '/home/users/qd201969/uor_track/fig'

lonl=0  #0  #
lonr=150#360#
lats=15  #
latn=70 #

def main_run():

def calc_monthly_max_mean(nl,nint):
    # nint, intensity columns in fftadd file
    # nint=9,max10mwind; nint=10,average precip
    ds = xr.open_dataset("%s/grid_%.1f.nc"%(fileout,radiu))
    grid = ds['grid']
    dim = grid.shape
    mean = np.zeros([12,dim[0],dim[0]],dtype=float)
    maxv = np.zeros([12,dim[0],dim[0]],dtype=float)
    numb = np.zeros([12,dim[0],dim[0]],dtype=float)

    ctime,clat,clon,cinte = composite_time("%s/%s_%d_1980-2020%s"%(
        path,prefix,nl,suffix),lats,latn,lonl,lonr,nint)
    for ct,clo,cla,cit in zip(ctime,clon,clat,cinte):




def gene_grid():
    lat = np.arange(lats,latn+0.1,radiu)
    lon = np.arange(lonl,lonr+0.1,radiu)
    var = np.zeros([len(lat),len(lon)],dtype=float)
    da = xr.DataArray(var, coords=[lat,lon], dims=["lat", "lon"])
    ds = da.to_dataset(name="grid")
    ds.to_netcdf("%s/grid_%.1f.nc"%(fileout,radiu),"w")

if __name__=='__main__':
    main_run()
