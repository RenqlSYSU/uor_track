import xarray as xr
import pandas as pd
import numpy as np
import subprocess
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import colors
import cartopy.crs as ccrs
from wrf import (to_np, getvar, interplevel,ll_to_xy,xy_to_ll)

def read_wrf_height(filname,hgt,lats,lonl): 
    ncfile = Dataset(filname)
    #x_y = to_np(ll_to_xy(ncfile,[llats,llatn],[llonl,llonr])) # return x and y
    u = getvar(ncfile, 'U10')
    v = getvar(ncfile, 'V10')
    term = u
    term.data = np.power((np.power(u.data,2)+np.power(v.data,2)),0.5)
    del u, v
    #term= getvar(ncfile, "wspd_wdir")[0]
    #hg  = getvar(ncfile, "height") #Model Height for Mass Grid (AGL)
    x_y = to_np(ll_to_xy(ncfile,lats,lonl)) # return x and y
    #var = interplevel(term, hg, hgt)[x_y[1],x_y[0]]
    var = term.sel(south_north=x_y[1],west_east=x_y[0])
    #del hg, x_y, term
    print(var)
    return var

# constants
BIGFONT=22
MIDFONT=18
SMFONT=14

lat = 22.1817
lon = 113.5600
hgt = 40 # unit m
lats = 22
latn = 23
lonl = 113
lonr = 114

# define date of start time and forcast time
stime  = ['2021','09','15'] # year, month, date, hour
ftime  = pd.date_range(start='2021-09-16 00',end='2021-09-16 23',freq='1H',closed=None)
ftime0 = pd.date_range(start='2021-09-16 00',end='2021-09-25 23',freq='1H',closed=None)

case = ['Coupled','Uncoupled',"Obs"]
path = ['/home/lzhenn/cooperate/data/Njord/njord/',
        '/home/dataop/data/nmodel/wrf_fc/',
        "/home/lzhenn/cooperate/data/case_study/210916/A_WIND_221817_1135600.csv"]
domin = ['d03','d04'] #njord, path
figdir= "/home/lzhenn/cooperate/fig/ts.jpg"
title = "40m wind, lat=22.1817, lon=113.5600 [MC_PS]"

var = np.empty((3,len(ftime)), dtype = float)
var[2,:] = np.loadtxt(path[2], skiprows=4, delimiter=",", usecols=2)
for nt in range(len(ftime)):
    for nc in range(len(case)-1):
        if case[nc] == case[0]: 
            filedir = path[nc]+''.join(stime)+\
                    '/wrfout_'+domin[nc]+'_'+ftime[nt].strftime("%Y-%m-%d_%H:00:00")
            filedir0= path[nc]+''.join(stime)+\
                    '/wrfout_'+domin[nc]+'_'+ftime0[nt].strftime("%Y-%m-%d_%H:00:00")
        else:
            filedir = path[nc]+stime[0]+'/'+''.join(stime[0:2])+'/'+''.join(stime)+\
                    '12/wrfout_'+domin[nc]+'_'+ftime[nt].strftime("%Y-%m-%d_%H:00:00")
            filedir0= path[nc]+stime[0]+'/'+''.join(stime[0:2])+'/'+''.join(stime)+\
                    '12/wrfout_'+domin[nc]+'_'+ftime0[nt].strftime("%Y-%m-%d_%H:00:00")

        var[nc,nt] = read_wrf_height(filedir,hgt,lat,lon) 
print(var)

fig = plt.figure(figsize=(12,9),dpi=300)
axe = fig.add_axes([0.05, 0.05, 0.9, 0.65])
axe.set_title(title,fontsize=MIDFONT) 
for nc in range(len(var)):
    axe.plot(ftime, var[nc,:], label=case[nc])
axe.set_ylabel("wind speed (m/s)",fontsize=MIDFONT)  # Add an x-label to the axes.
axe.set_xlabel("time",fontsize=MIDFONT)  # Add a y-label to the axes.
axe.tick_params(axis='both', which='major', labelsize=SMFONT)
#axe.set_xticklabels(fontsize=SMFONT)
#axe.set_yticklabels(fontsize=SMFONT)
axe.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d_%H:00"))
plt.setp(axe.get_xticklabels(), rotation=30, ha="right")
axe.legend(fontsize=MIDFONT)  # Add a legend.
fig.savefig(figdir,bbox_inches="tight",pad_inches=0.1)
#fig.savefig("%s%s.jpg"%(figdir,title),bbox_inches="tight",pad_inches=0.1)

