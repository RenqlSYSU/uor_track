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

indir = ['/home/lzhenn/array74/data/ssp_rainfall/hist_2010s/',
         '/home/metctm1/array/data/imerg/mon/3B-MO.MS.MRG.3IMERG.']
case = ['hist_2010s','IMERG']
years = range(2011,2021)
maskwater = True

ncfile = Dataset('/home/lzhenn/cmip6-wrf-arch/hist_2010s/2011/'+
    'analysis/wrfout_%s_2011-01-02_04:00:00'%dom)
landmask = wrf.getvar(ncfile,'LANDMASK').data
print(landmask)
#LAND MASK (1 FOR LAND, 0 FOR WATER)

lat_sp = 1.0 #5.0
lon_sp = 1.0 #10.0
    
def main_run():
    dom = ['d02','d03','d04']
    perid = ['Apr-Sep','AMJ','JAS']


def draw_precip_daily():
    fig = plt.figure(figsize=(12,9),dpi=100)
    ax = fig.add_axes([0.1, 0.1, 0.9, 0.65])
    ax.set_title("%s Daily cycle"%dom,fontsize=title_font,fontdict=font)
   
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

if __name__=='__main__':
    main_run()
