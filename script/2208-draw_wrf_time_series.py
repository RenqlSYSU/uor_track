import xarray as xr
import pandas as pd
import numpy as np
import gc #garbage collector
import subprocess, os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import colors
from matplotlib import cm
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.io.shapereader import Reader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmaps, wrf
from netCDF4 import Dataset
from datetime import datetime

title_font=14
label_font=10
plt.rcParams["font.weight"] = "bold"
font = {'family': 'sans-serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black',
        }

figdir = "/home/lzhenn/cooperate/fig"
case = ['WRFONLY','WRFROMS']
indir = ['/home/metctm1/array_hq133/data/s2s/wrfonly/clim/',
         '/home/metctm1/array_hq133/data/s2s/wrfroms/clim/aegir_']
lats = 30  #15  #36  #
latn = 31  #19  #37  #
lonl = 120 #110 #100 #
lonr = 121 #118 #101 #

def main_run():
    draw_wrf_ts('SST')
    #draw_wrf_ts_1grid('TSK')
    #draw_wrf_ts_1grid('T2')

def draw_wrf_ts_1grid(varname):
    wrf1=Dataset('%s/2016070100/wrfout_d01_2016-07-23_00:00:00'%indir[0])
    lake = wrf.getvar(wrf1,'LAKEMASK').data #1 FOR LAKE, 0 FOR NON-LAKE
    lat = wrf.getvar(wrf1,'XLAT').data
    lon = wrf.getvar(wrf1,'XLONG').data
    ind = np.argwhere(np.array([lat<=latn,lat>=lats,lon<=lonr,lon>=lonl,
        lake==1]).all(axis=0))
    print(ind.shape)
    ni = int(len(ind)/2)
    print(ni)
    ilat = lat[ind[ni,0],ind[ni,1]]
    ilon = lon[ind[ni,0],ind[ni,1]]
    del wrf1, lat, lon, lake 
    gc.collect()

    var = [] 
    for nc in range(len(indir)):
        fn_stream = subprocess.check_output(
            'ls %s2011070100/wrfout_d01_2011-07-*'%(
            indir[nc]), shell=True).decode('utf-8')
        fn_list = fn_stream.split()
        
        wrf_list=[Dataset(itm) for itm in fn_list]
        print('total file %d; opened file %d'%(
            len(fn_list),len(wrf_list)))

        var.append(wrf.getvar(wrf_list, varname, timeidx=wrf.ALL_TIMES, 
            method='cat').data[:,ind[ni,0],ind[ni,1]]-273.15)
        print(var)
        
        if nc==1:
            time = wrf.getvar(wrf_list, 'Times', timeidx=wrf.ALL_TIMES, 
                method='cat').data
            ftime = pd.to_datetime(np.datetime_as_string(time,
                unit='h'),format='%Y-%m-%dT%H')
            del time
        del wrf_list
        gc.collect()

    draw_ts('%s (%.2fN, %.2fE)'%(varname,ilat,ilon),
        ftime,var,'%s/%s_ts.png'%(figdir,varname))
    
def draw_wrf_ts(varname):
    wrf1=Dataset('%s/2016070100/wrfout_d01_2016-07-23_00:00:00'%indir[0])
    #x_y = wrf.ll_to_xy(wrf1,[lats,latn],[lonl,lonr]).data
    #print(x_y)
    lat = wrf.getvar(wrf1,'XLAT').data
    lon = wrf.getvar(wrf1,'XLONG').data
    ind = np.argwhere(np.array([lat<=latn,lat>=lats,lon<=lonr,lon>=lonl]
        ).all(axis=0))
    del wrf1, lat, lon 
    gc.collect()
    print(ind.shape)

    var = [] 
    for nc in range(len(indir)):
        fn_stream = subprocess.check_output(
            'ls %s2016070100/wrfout_d01_2016-08-*'%(
            indir[nc]), shell=True).decode('utf-8')
        fn_list = fn_stream.split()
        
        wrf_list=[Dataset(itm) for itm in fn_list]
        print('total file %d; opened file %d'%(
            len(fn_list),len(wrf_list)))

        #term = wrf.getvar(wrf_list, varname, timeidx=wrf.ALL_TIMES, 
        #    method='cat')[:,x_y[1,0]:x_y[1,1],x_y[0,0]:x_y[0,1]]
        term = wrf.getvar(wrf_list, varname, timeidx=wrf.ALL_TIMES, 
            method='cat').data[:,ind[:,0],ind[:,1]]
        print(term.shape)
        #var.append(np.nanmean(term.data,(1,2))-273.15)
        var.append(np.nanmean(term,1)-273.15)
        print(var)

    time = wrf.getvar(wrf_list, 'Times', timeidx=wrf.ALL_TIMES, 
        method='cat').data
    ftime = pd.to_datetime(np.datetime_as_string(time,
        unit='h'),format='%Y-%m-%dT%H')
    del wrf_list, time
    gc.collect()
    draw_ts(varname,ftime,var,'%s/%s_ts.png'%(figdir,varname))
   
def draw_ts(figtile,ftime,var,outdir):
    fig = plt.figure(figsize=(12,4),dpi=150) #width,height
    axe = fig.subplots(1,1)
    axe.set_title('%s'%figtile,fontsize=title_font,fontdict=font)
    frmt = ['ro-','bo-']
    for nr in range(len(var)):
        axe.plot(ftime, var[nr], frmt[nr], linewidth=2, markersize=2)
    
    axe.set_ylabel('',fontsize=label_font,fontdict=font)
    axe.set_xlabel('',fontsize=label_font,fontdict=font)
    axe.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d %H'))
    plt.setp(axe.get_xticklabels(), rotation=15, ha="right")
    axe.legend(case)
    
    plt.savefig(outdir,bbox_inches='tight',pad_inches=0.01)

if __name__=='__main__':
    main_run()

