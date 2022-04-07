import xarray as xr
import pandas as pd
import numpy as np
import os, subprocess, wrf
from multiprocessing import Pool
from netCDF4 import Dataset
from datetime import datetime, timedelta
import matplotlib.pyplot as plt

title_font=18
label_font=14
plt.rcParams["font.weight"] = "bold"
font = {'family': 'serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black',
        }

indir = '/home/lzhenn/cmip6-wrf-arch/'
otdir = '/home/lzhenn/array74/data/ssp_rainfall/'
case = ['hist_2010s','ssp245_2040s_wholeyear','ssp245_2090s_wholeyear']
years= [range(2011,2021),range(2040,2050),range(2090,2100)]
dom = 'd02'

def main_run():
    #bins = [0.1,0.5,1,2,3,5,8,12,16,20,25,30,40,50,70,100,150] # hourly preci
    #bins = [0.1,0.5,1,3,5,8,12,16,20,25,30,40,50,70,100,150,200] # 3h preci
    bins = [0.1,0.5,1,3,5,10,15,20,30,40,60,80,100,120,150,200,250] #24h accumulated preci
    for i in range(len(case)):
        calc_monthly_preci(i,case[i],years[i])
    #    draw_precip_pdf(case[i],years[i],24,bins)
    '''
    ntask = len(case)
    process_pool = Pool(processes=ntask)
    results=[]
    for i in range(0,ntask):
        result=process_pool.apply_async(calc_monthly_preci,
            args=(i,case[i],years[i]))
        results.append(result)
    print(results)
    print(results[0].get())
    process_pool.close()
    process_pool.join()
    print(results[0].get())
    '''

def calc_monthly_preci(i,case,years):
    for year in years:
        fn_stream = subprocess.check_output(
            'ls %s%s/%d/analysis/wrfout_%s_%d-*'%(
            indir,case,year,dom,year), shell=True).decode('utf-8')
        fn_list = np.array(fn_stream.split())
        iyear  = np.array([int(a.rsplit('_')[-2].split('-')[0]) for a in fn_list])
        imonth = np.array([int(a.rsplit('_')[-2].split('-')[1]) for a in fn_list])

        for nm in range(1,13):
            datafile = "%s%s/wrfout_%s.precc.%d%02d.nc"%(otdir,case,dom,year,nm)
            if not os.path.exists(datafile):
                indx = np.argwhere(np.array([iyear==year,imonth==nm]).all(axis=0))
                calc_hourly_preci(i,case,year,nm,datafile,fn_list[indx[:,0]])

def calc_hourly_preci(i,case,year,nm,datafile,fn_list):
    wrf_list=[Dataset(itm) for itm in fn_list]
    ntime = len(wrf_list)
    print("task[%d] %s %d-%2d time frames: %d"%(i,case,year,nm,ntime))

    rainnc = wrf.getvar(wrf_list,'RAINNC',timeidx=wrf.ALL_TIMES, 
        method='cat').data
    rainc = wrf.getvar(wrf_list,'RAINC',timeidx=wrf.ALL_TIMES, 
        method='cat').data
    var = rainnc + rainc
    var[1:ntime,:,:] = var[1:ntime,:,:] - var[0:ntime-1,:,:]
    
    time = wrf.getvar(wrf_list,'Times',timeidx=wrf.ALL_TIMES, 
        method='cat').data
    time = pd.to_datetime(np.datetime_as_string(
        time,unit='h'),format='%Y-%m-%dT%H')
    lat = wrf.getvar(wrf_list[0],'XLAT').data
    lon = wrf.getvar(wrf_list[0],'XLONG').data
   
    if nm > 1:
        a = (time[0]-timedelta(hours=1)).strftime("%Y-%m-%d_%H:00:00")
        print(fn_list[0])
        print('%s%s/%d/analysis/wrfout_%s_%s'%(indir,case,year,dom,a))
        ncfile = Dataset('%s%s/%d/analysis/wrfout_%s_%s'%(
            indir,case,year,dom,a))
        b = wrf.getvar(ncfile,"RAINC").values
        c = wrf.getvar(ncfile,"RAINNC").values
        var[0,:,:] = var[0,:,:] - (b+c)
    
    da = xr.DataArray(var,coords=[time,lat[:,0],lon[0,:]],dims=['time','lat','lon'])
    da.attrs['units']='mm/h'
    ds = da.to_dataset(name="temp")
    ds.to_netcdf(datafile,"w")

def draw_precip_pdf(case,years,hr,bins):
    datafile = "%s%s/wrfout_%s.precc.*.nc"%(otdir,case,dom)
    ds = xr.open_mfdataset(datafile,concat_dim='time', combine='nested')
    var = ds['temp'].data
    dim = np.array(var.shape)
    print(dim)
    
    if hr>1:
        term = [var[i:i+hr,:,:].sum(axis=0) for i in range(0,dim[0],hr)]
        dim[0]=dim[0]/hr
    else:
        term = var
    var = np.array(term).reshape(np.prod(dim))
    print(dim)
    #numb, bins = np.histogram(var,20)
    #numb = numb*100/np.prod(dim)
    
    fig = plt.figure(figsize=(12,9),dpi=100)
    ax = fig.add_axes([0.1, 0.1, 0.9, 0.65])
    n,binss,patches = ax.hist(var,bins,density=True,log=True)
    print(n)
    
    ax.set_ylim(10e-6,1)
    ax.grid(True, which="both", color='grey', linestyle='--', linewidth=1)
    ax.set_xlabel('precip (mm/%dh)'%hr,fontsize=label_font,fontdict=font)
    ax.set_ylabel("probability density",fontsize=label_font,fontdict=font)
    ax.set_title("%s"%case,fontsize=title_font,fontdict=font)
    
    fig.savefig("/home/lzhenn/cooperate/fig/cmip6_%s_precip_%dh"%(
        case,hr),bbox_inches='tight',pad_inches=0.1)
    
if __name__=='__main__':
    main_run()
