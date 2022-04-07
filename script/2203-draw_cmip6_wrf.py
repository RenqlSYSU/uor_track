import xarray as xr
import pandas as pd
import numpy as np
import os, subprocess, wrf
from multiprocessing import Pool
from netCDF4 import Dataset
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns

title_font=18
label_font=14
plt.rcParams["font.weight"] = "bold"
font = {'family': 'serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black',
        }

otdir = '/home/lzhenn/array74/data/ssp_rainfall/'
cases = ['hist_2010s','ssp245_2040s_wholeyear','ssp245_2090s_wholeyear']
label = ['2010s','2040s','2090s']
years = [range(2011,2021),range(2040,2050),range(2090,2100)]
dom = 'd04'
maskwater = True

ncfile = Dataset('/home/lzhenn/cmip6-wrf-arch/hist_2010s/2011/'+
    'analysis/wrfout_%s_2011-01-02_04:00:00'%dom)
landmask = wrf.getvar(ncfile,'LANDMASK').data
print(landmask)
#LAND MASK (1 FOR LAND, 0 FOR WATER)

def main_run():
    #bins = [0.1,0.5,1,2,3,5,8,12,16,20,25,30,40,50,70,100,150] # hourly preci
    #bins = [0.1,0.5,1,3,5,8,12,16,20,25,30,40,50,70,100,150,200] # 3h preci
    bins = [0.1,0.5,1,3,5,10,15,20,30,40,60,80,100,120,150,200,250] #24h accumulated preci
    #draw_precip_pdf_3case(24, bins)
    draw_precip_annual()
    #draw_precip_daily()
    #draw_precip_monthly_box(24)
    #draw_precip_monthly_box(3)
    #draw_precip_monthly_box(1)

def draw_precip_monthly_box(hr):
    fig = plt.figure(figsize=(12,9),dpi=100)
    ax = fig.add_axes([0.1, 0.1, 0.9, 0.65])
    ax.set_title("%s Precip (mm/%dhr)"%(hr,dom),fontsize=title_font,fontdict=font)
   
    dicts = {'case':[],'month':[],'precip':[]}
    for case,labl in zip(cases,label):
        for nm in range(1,13):
            datafile = "%s%s/wrfout_%s.precc.*%02d.nc"%(otdir,case,dom,nm)
            ds = xr.open_mfdataset(datafile,concat_dim='time', combine='nested')
            var = ds['temp'].data
            dim = np.array(var.shape)
            if hr>1:
                term = [var[i:i+hr,:,:].sum(axis=0) for i in range(0,dim[0],hr)]
                dim[0]=dim[0]/hr
            else:
                term = var
            var = list(np.array(term).reshape(np.prod(dim)))
            dicts['case']  = dicts['case']+[case]*len(var)
            dicts['month'] = dicts['month']+[nm]*len(var)
            dicts['precip']= dicts['precip']+var
    df = pd.DataFrame(dicts)
    print(df)

    sns.boxplot(x='month', y='precip', hue="case", data=df, palette="Set1")
    ax.set_ylim([0,250])
    fig.savefig("/home/lzhenn/cooperate/fig/cmip6_precip_monthly_box_%dh"%hr
        ,bbox_inches='tight',pad_inches=0.1)

def draw_precip_daily():
    fig = plt.figure(figsize=(12,9),dpi=100)
    ax = fig.add_axes([0.1, 0.1, 0.9, 0.65])
    ax.set_ylabel('precip (mm/h)',fontsize=label_font,fontdict=font)
    ax.set_xlabel("hour",fontsize=label_font,fontdict=font)
    ax.set_title("%s Daily cycle"%dom,fontsize=title_font,fontdict=font)
   
    for season,shape in zip(['AMJ','JAS'],['-','-.']):
        for case,labl,colr in zip(cases,label,['black','red','blue']):
            if season == 'AMJ':
                datafile = "%s%s/wrfout_%s.precc.*0[456].nc"%(otdir,case,dom)
            else:
                datafile = "%s%s/wrfout_%s.precc.*0[789].nc"%(otdir,case,dom)
            ds = xr.open_mfdataset(datafile,concat_dim='time', combine='nested')
            var = ds['temp']
            if maskwater:
                mask = np.tile(landmask,(len(var),1,1))==0
                var.data = np.ma.array(var.data,mask=mask)
            var = var.mean(('lat','lon'))
            var = var.groupby(var.time.dt.hour).mean('time')
            print(var)
            ax.plot(var.hour,var.data,linestyle=shape,linewidth=2,
                color=colr,label='%s-%s'%(labl,season))
    
    ax.tick_params(axis='both', which='major', labelsize=label_font)
    ax.grid(True, which="both", color='grey', linestyle='--', linewidth=1)
    ax.legend(loc='upper right') 
    fig.savefig("/home/lzhenn/cooperate/fig/cmip6_precip_daily"
        ,bbox_inches='tight',pad_inches=0.1)

def draw_precip_annual():
    dtime = pd.date_range(start='2019-01-01',end='2019-12-31',freq='5d',closed=None)

    fig = plt.figure(figsize=(12,9),dpi=100)
    ax = fig.add_axes([0.1, 0.1, 0.9, 0.65])
    ax.set_ylabel('precip (mm/5d)',fontsize=label_font,fontdict=font)
    ax.set_title("%s Annual cycle"%dom,fontsize=title_font,fontdict=font)
    
    for case,labl in zip(cases,label):
        datafile = "%s%s/wrfout_%s.precc.*.nc"%(otdir,case,dom)
        ds = xr.open_mfdataset(datafile,concat_dim='time',combine='nested')
        var = ds['temp']
        if maskwater:
            mask = np.tile(landmask,(len(var),1,1))==0
            var.data = np.ma.array(var.data,mask=mask)
        var = var.mean(('lat','lon'))
        var = var.resample(time="24H").sum()
        var = var.groupby(var.time.dt.dayofyear).mean('time')
        term = [var[i:i+5].sum(axis=0) for i in range(0,365,5)]
        term = np.array(term)
        print(var)
        ax.plot(dtime,term,linewidth=2,label=labl)
    
    ax.tick_params(axis='both', which='major', labelsize=label_font)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    ax.grid(True, which="both", color='grey', linestyle='--', linewidth=1)
    ax.legend(loc='upper right') 
    fig.savefig("/home/lzhenn/cooperate/fig/cmip6_precip_annual"
        ,bbox_inches='tight',pad_inches=0.1)

def draw_precip_pdf_3case(hr,bins):
    fig = plt.figure(figsize=(12,9),dpi=100)
    ax = fig.add_axes([0.1, 0.1, 0.9, 0.65])
    ax.set_xlabel('precip (mm/%dh)'%hr,fontsize=label_font,fontdict=font)
    ax.set_ylabel("probability density",fontsize=label_font,fontdict=font)
    ax.set_title("%s Precip pdf"%dom,fontsize=title_font,fontdict=font)
    
    for case,labl in zip(cases,label):
        datafile = "%s%s/wrfout_%s.precc.*.nc"%(otdir,case,dom)
        ds = xr.open_mfdataset(datafile,concat_dim='time', combine='nested')
        var = ds['temp'].data
        dim = np.array(var.shape)
        print('%s dimension '%case,dim)
        
        if hr>1:
            term = [var[i:i+hr,:,:].sum(axis=0) for i in range(0,dim[0],hr)]
            dim[0]=dim[0]/hr
        else:
            term = var
        
        var = np.array(term).reshape(np.prod(dim))
        print('%s dimension '%case,dim)
        numb, bins = np.histogram(var,bins)
        numb = numb/np.prod(dim)
        ax.plot(bins[1:],numb,linewidth=2,label=labl)
    
    ax.set_yscale('log')
    ax.set_xscale('log')
    #ax.set_ylim(10e-6,1)
    ax.grid(True, which="both", color='grey', linestyle='--', linewidth=1)
    ax.legend(loc='upper right') 
    
    fig.savefig("/home/lzhenn/cooperate/fig/cmip6_precip_%dh_pdf"%(
        hr),bbox_inches='tight',pad_inches=0.1)

if __name__=='__main__':
    main_run()
