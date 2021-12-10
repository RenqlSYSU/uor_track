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
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import colors
import cartopy.crs as ccrs
from wrf import (to_np, getvar, interplevel,ll_to_xy,xy_to_ll,ALL_TIMES)

def read_wrf_height_rh(ncfile,hgt,lats,lonl): 
    term= getvar(ncfile, "rh", timeidx=ALL_TIMES, method="cat")
    hg  = getvar(ncfile, "height", timeidx=ALL_TIMES, method="cat") #Model Height for Mass Grid (AGL)
    x_y = to_np(ll_to_xy(ncfile,lats,lonl)) # return x and y
    var = interplevel(term, hg, hgt)[:,x_y[1],x_y[0]]
    print(term)
    del hg, x_y, term
    return var

def read_wrf_height_temp(ncfile,hgt,lats,lonl): 
    term= getvar(ncfile, "temp", units="degC", timeidx=ALL_TIMES, method="cat")
    hg  = getvar(ncfile, "height", timeidx=ALL_TIMES, method="cat") #Model Height for Mass Grid (AGL)
    x_y = to_np(ll_to_xy(ncfile,lats,lonl)) # return x and y
    var = interplevel(term, hg, hgt)[:,x_y[1],x_y[0]]
    print(term)
    del hg, x_y, term
    return var

def read_wrf_height_wnd(ncfile,hgt,lats,lonl,index=0): 
    #x_y = to_np(ll_to_xy(ncfile,[llats,llatn],[llonl,llonr])) # return x and y
    #u = getvar(ncfile, 'U10')
    #v = getvar(ncfile, 'V10')
    #term = u
    #term.data = np.power((np.power(u.data,2)+np.power(v.data,2)),0.5)
    #del u, v
    #var = term.sel(south_north=x_y[1],west_east=x_y[0])
    term= getvar(ncfile, "wspd_wdir", timeidx=ALL_TIMES, method="cat")[index] # index = 0 speed, 1 direction 
    hg  = getvar(ncfile, "height", timeidx=ALL_TIMES, method="cat") #Model Height for Mass Grid (AGL)
    x_y = to_np(ll_to_xy(ncfile,lats,lonl)) # return x and y
    var = interplevel(term, hg, hgt)[:,x_y[1],x_y[0]]
    del hg, x_y, term
    return var

# constants
BIGFONT=22
MIDFONT=20
SMFONT=18

lat = 20.961667 
lon = 111.608889
hgt = [110,110,50,50,30,30,20,20] # unit m
drawvar = ['WS','WD','WS','WD','WS','WD','Temp','RH']
unit = ['m/s','°','m/s','°','m/s','°','°C','%']
#unit = ['m/s','Ã‚Â°','m/s','Ã‚Â°','m/s','Ã‚Â°','Ã‚Â°C','%']

# define date of start time and forcast time
stime  = ['2019','11','30'] # year, month, date, hour
ftime  = pd.date_range(start='2019-12-01 00',end='2019-12-31 00',freq='1H',closed=None)
ftime0 = pd.date_range(start='2021-09-16 00',end='2021-09-25 23',freq='1H',closed=None) # use for rainfall

case = ["Obs","Coupled"]#,'Uncoupled'
path = ["/home/lzhenn/array74/data/yangjiang-windfarm/1hour.txt",
        '/home/lzhenn/cooperate/data/case_study/yangjiang-windfarm/'+''.join(stime)+'00/']
        #'/home/lzhenn/cooperate/data/case_study/yangjiang-windfarm/'+''.join(stime)+'00_wrf/',
domin = ['d02','d02'] #njord, pathn2

for nv in range(0,8):
    varname = "%s_%dm_Avg [%s]"%(drawvar[nv],hgt[nv],unit[nv])
    figdir= "/home/lzhenn/cooperate/fig/ts_yangjiang1_%s.jpg"%(varname.rsplit("_",1)[0])
    title = "%s, %fN, %fE"%(varname,lat,lon)

    var = np.zeros((len(case),len(ftime)), dtype = float)
    df = pd.read_csv(path[0],index_col=0, usecols=['Date/Time',varname])
    var[0,:] = df[df.index.isin(ftime.strftime('%Y-%m-%d %H:00'))].to_numpy()[:,0]
    #df = pd.read_excel(path[0], index_col=0, usecols=['日期',varname])
    #var[0,:] = df[df.index.strftime('%Y-%m-%d %H').isin(ftime.strftime('%Y-%m-%d %H'))].to_numpy()[:,0]
    del df

    for nc in range(1,len(case),1):
        for nt in range(0,len(ftime)-1,24):
            #filedir = path[nc]+'wrfout_'+domin[nc]+'_'+ftime[nt].strftime("%Y-%m-%d_%H:00:00")
            filedir = path[nc]+'wrfout_'+domin[nc]+'_'+ftime[nt].strftime("%Y-%m-%d")+"*"
            fn_stream = subprocess.check_output('ls '+filedir, shell=True).decode('utf-8')
            fn_list   = fn_stream.split()
            fillist = []
            for itm in fn_list:
                print(itm)
                fillist.append(Dataset(itm))
            if drawvar[nv] == "WD":
                var[nc,nt:(nt+24)] = read_wrf_height_wnd(fillist,hgt[nv],lat,lon,index=1) 
            elif drawvar[nv] == "WS":
                var[nc,nt:(nt+24)] = read_wrf_height_wnd(fillist,hgt[nv],lat,lon,index=0) 
            elif drawvar[nv] == "Temp":
                #var[nc,nt:(nt+24)] = read_wrf_height_temp(fillist,hgt[nv],lat,lon) 
                var[nc,nt:(nt+24)] = read_wrf_height_temp(fillist,40,lat,lon) 
            elif drawvar[nv] == "RH":
                #var[nc,nt:(nt+24)] = read_wrf_height_rh(fillist,hgt[nv],lat,lon) 
                var[nc,nt:(nt+24)] = read_wrf_height_rh(fillist,40,lat,lon) 

    fig = plt.figure(figsize=(12,9),dpi=300)
    axe = fig.add_axes([0.05, 0.05, 0.9, 0.45])
    axe.set_title(title,fontsize=MIDFONT) 
    for nc in range(len(case)):
        if nc == 0:
            axe.plot(ftime, var[nc,:], label=case[nc])
        else:
            rmse = np.sqrt(np.power(np.nan_to_num((var[nc,:] - var[0,:]),copy=False,nan=0),2).mean())
            mae = np.nan_to_num((var[nc,:] - var[0,:]),copy=False,nan=0).mean()
            axe.plot(ftime, var[nc,:], label="%s (RMSE: %.2f, MAE: %.2f)"%(case[nc],rmse,mae))
    #axe.set_ylabel("wind speed (m/s)",fontsize=MIDFONT)  # Add an x-label to the axes.
    axe.set_xlabel("time",fontsize=MIDFONT)  # Add a y-label to the axes.
    axe.tick_params(axis='both', which='major', labelsize=SMFONT)
    axe.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d_%H:00"))
    plt.setp(axe.get_xticklabels(), rotation=30, ha="right")
    axe.legend(fontsize=MIDFONT)  # Add a legend.
    fig.savefig(figdir,bbox_inches="tight",pad_inches=0.1)

