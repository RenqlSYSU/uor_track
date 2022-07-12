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
cases = ['hist_2010s','ssp245_2040s_wholeyear','ssp245_2090s_wholeyear','IMERG']
label = ['2010s','2040s','2090s','IMERG']
years = [range(2011,2021),range(2040,2050),range(2090,2100)]
maskwater = True

dom = 'd04'
ncfile = Dataset('/home/lzhenn/cmip6-wrf-arch/hist_2010s/2011/'+
    'analysis/wrfout_%s_2011-01-02_04:00:00'%dom)
landmask = wrf.getvar(ncfile,'LANDMASK').data
flat = wrf.getvar(ncfile,'XLAT').data[:,0]
flon = wrf.getvar(ncfile,'XLONG').data[0,:]
print(landmask)
print(flat)
print(flon)
#LAND MASK (1 FOR LAND, 0 FOR WATER)

def main_run():
    #bins = [0.1,0.5,1,2,3,5,8,12,16,20,25,30,40,50,70,100,150] # hourly preci
    #bins = [0.1,0.5,1,3,5,8,12,16,20,25,30,40,50,70,100,150,200] # 3h preci
    bins = [0.1,0.5,1,3,5,10,15,20,30,40,60,80,100,120,150,200,250] #24h accumulated preci
    #draw_precip_pdf_3case(24, bins, dom)
    draw_precip_annual(dom)
    draw_precip_daily(dom)
    #draw_precip_monthly_box(24,dom)

def draw_precip_monthly_box(hr,dom):
    fig = plt.figure(figsize=(12,9),dpi=100)
    ax = fig.add_axes([0.1, 0.1, 0.9, 0.65])
    ax.set_title("%s Precip (mm/%dhr) NoOcean=%s"%(dom,hr,maskwater),fontsize=title_font,fontdict=font)
   
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
            
            if maskwater:
                mask = np.tile(landmask,(len(term),1,1))==0
                var = list(np.ma.array(np.array(term),mask=mask).max(axis=(0)).compressed())
                #var = list(np.ma.array(np.array(term),mask=mask).mean(axis=(1,2)))
                print('dimension %d'%len(var))
                #var = list(np.ma.array(np.array(term),mask=mask).compressed())
            else:
                var = list(np.array(term).reshape(np.prod(dim)))
            
            dicts['case']  = dicts['case']+[case]*len(var)
            dicts['month'] = dicts['month']+[nm]*len(var)
            dicts['precip']= dicts['precip']+var
    df = pd.DataFrame(dicts)
    print(df)

    sns.boxplot(x='month', y='precip', hue="case", data=df, palette="Set1")
    #ax.set_ylim([0,250])
    fig.savefig("/home/lzhenn/cooperate/fig/cmip6_precip_monthly_box_%dh_%s"%(hr,dom)
        ,bbox_inches='tight',pad_inches=0.1)

def draw_precip_daily(dom):
    fig = plt.figure(figsize=(12,9),dpi=100)
    ax = fig.add_axes([0.1, 0.1, 0.9, 0.65])
    ax.set_ylabel('precip (mm/h)',fontsize=label_font,fontdict=font)
    ax.set_xlabel("hour",fontsize=label_font,fontdict=font)
    ax.set_title("%s Daily cycle NoOcean=%s"%(dom,maskwater),fontsize=title_font,fontdict=font)
   
    for season,shape in zip(['AMJ','JAS'],['-','-.']):
        for case,labl,colr in zip(cases,label,['black','tab:blue','tab:orange','red']):
            if case == 'IMERG':
                path = '/home/metctm1/array/data/imerg/3hr'
                datafile = '%s/2011/3B-HHR.MS.MRG.3IMERG.20110423-*.nc4'%(path)
                ds = xr.open_mfdataset(datafile,concat_dim='time', combine='nested')
                lat = ds.lat
                lon = ds.lon
                ilon = lon[(lon>=flon[0]) & (lon<=flon[-1])]
                ilat = lat[(lat>=flat[0]) & (lat<=flat[-1])]
                var = ds['precipitationCal'].sel(lat=ilat,lon=ilon)
                var.data = np.zeros(var.shape)
                for year in years[0]:
                    if season == 'AMJ':
                        datafile = '%s/%d/3B-HHR.MS.MRG.3IMERG.%d0[456]*.nc4'%(path,year,year)
                    else:
                        datafile = '%s/%d/3B-HHR.MS.MRG.3IMERG.%d0[789]*.nc4'%(path,year,year)
                    ds = xr.open_mfdataset(datafile,concat_dim='time', combine='nested')
                    term = ds['precipitationCal'].sel(lat=ilat,lon=ilon)
                    var.data = var.data + term.groupby(term.time.dt.hour).mean('time')
                var.data = var.data/len(years[0])
                var = var.transpose('time','lat','lon')
                print('mask water')
                var = var.interp(lat=flat,lon=flon)
                if maskwater:
                    mask = np.tile(landmask,(len(var),1,1))==0
                    var.data = np.ma.array(var.data,mask=mask)
            else:
                if season == 'AMJ':
                    datafile = "%s%s/wrfout_%s.precc.*0[456].nc"%(otdir,case,dom)
                else:
                    datafile = "%s%s/wrfout_%s.precc.*0[789].nc"%(otdir,case,dom)
                ds = xr.open_mfdataset(datafile,concat_dim='time', combine='nested')
                var = ds['temp']
                var = var.groupby(var.time.dt.hour).mean('time')
                if maskwater:
                    mask = np.tile(landmask,(len(var),1,1))==0
                    var.data = np.ma.array(var.data,mask=mask)
            var = var.mean(('lat','lon'))
            print(var)
            ax.plot(range(0,24,1),var.data,linestyle=shape,linewidth=2,
                color=colr,label='%s-%s'%(labl,season))
    
    ax.tick_params(axis='both', which='major', labelsize=label_font)
    ax.grid(True, which="both", color='grey', linestyle='--', linewidth=1)
    ax.legend(loc='upper right') 
    fig.savefig("/home/lzhenn/cooperate/fig/cmip6_precip_daily_%s"%dom
        ,bbox_inches='tight',pad_inches=0.1)

def draw_precip_annual(dom):
    dtime = pd.date_range(start='2019-01-01',end='2019-12-31',freq='5d',closed=None)

    fig = plt.figure(figsize=(12,9),dpi=100)
    ax = fig.add_axes([0.1, 0.1, 0.9, 0.65])
    ax.set_ylabel('precip (mm/5d)',fontsize=label_font,fontdict=font)
    ax.set_title("%s Annual cycle NoOcean=%s"%(dom,maskwater),fontsize=title_font,fontdict=font)
    
    for case,labl,colr in zip(cases,label,['black','tab:blue','tab:orange','black']):
        if case == 'IMERG':
            datafile = "/home/metctm1/array/data/imerg/daily/3B-DAY.MS.MRG.3IMERG.*"
            ds = xr.open_mfdataset(datafile,concat_dim='time',combine='nested')
            lat = ds.lat
            lon = ds.lon
            ilon = lon[(lon>=flon[0]) & (lon<=flon[-1])]
            ilat = lat[(lat>=flat[0]) & (lat<=flat[-1])]
            var = ds['precipitationCal'].sel(lat=ilat,lon=ilon)
            var = var.transpose('time','lat','lon')
            var = var.groupby(var.time.dt.dayofyear).mean('time')
            if maskwater:
                print('mask water')
                var = var.interp(lat=flat,lon=flon)
                mask = np.tile(landmask,(len(var),1,1))==0
                var.data = np.ma.array(var.data,mask=mask)
            var = var.mean(('lat','lon'))
            term = [var[i:i+5].sum(axis=0) for i in range(0,365,5)]
            term = np.array(term)
            ax.plot(dtime,term,linewidth=2,label=labl,color=colr,linestyle='dotted')
        else:
            datafile = "%s%s/wrfout_%s.precc.*.nc"%(otdir,case,dom)
            ds = xr.open_mfdataset(datafile,concat_dim='time',combine='nested')
            var = ds['temp']
            var = var.resample(time="24H").sum()
            var = var.groupby(var.time.dt.dayofyear).mean('time')
            if maskwater:
                mask = np.tile(landmask,(len(var),1,1))==0
                var.data = np.ma.array(var.data,mask=mask)
            var = var.mean(('lat','lon'))
            term = [var[i:i+5].sum(axis=0) for i in range(0,365,5)]
            term = np.array(term)
            ax.plot(dtime,term,linewidth=2,label=labl,color=colr)
        
    ax.tick_params(axis='both', which='major', labelsize=label_font)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    ax.grid(True, which="both", color='grey', linestyle='--', linewidth=1)
    ax.legend(loc='upper right') 
    fig.savefig("/home/lzhenn/cooperate/fig/cmip6_precip_annual_%s"%dom
        ,bbox_inches='tight',pad_inches=0.1)

def draw_precip_pdf_3case(hr,bins,dom):
    fig = plt.figure(figsize=(12,9),dpi=100)
    ax = fig.add_axes([0.1, 0.1, 0.9, 0.65])
    ax.set_xlabel('precip (mm/%dh)'%hr,fontsize=label_font,fontdict=font)
    ax.set_ylabel("probability density",fontsize=label_font,fontdict=font)
    ax.set_title("%s Precip pdf NoOcean=%s"%(dom,maskwater),fontsize=title_font,fontdict=font)
    
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
        
        if maskwater:
            mask = np.tile(landmask,(len(term),1,1))==0
            var = np.ma.array(np.array(term),mask=mask).compressed()
        else:
            var = np.array(term).reshape(np.prod(dim))
        
        print('%s dimension '%case,dim)
        print('len(var) %d'%(len(var)))
        numb, bins = np.histogram(var,bins)
        numb = numb/len(var)
        ax.plot(bins[1:],numb,linewidth=2,label=labl)
    
    ax.set_yscale('log')
    ax.set_xscale('log')
    #ax.set_ylim(10e-6,1)
    ax.grid(True, which="both", color='grey', linestyle='--', linewidth=1)
    ax.legend(loc='upper right') 
    
    fig.savefig("/home/lzhenn/cooperate/fig/cmip6_precip_%dh_pdf_%s"%(
        hr,dom),bbox_inches='tight',pad_inches=0.1)

if __name__=='__main__':
    main_run()
