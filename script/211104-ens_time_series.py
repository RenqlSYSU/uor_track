"""
read Whole_HK_DTM_5m.asc, 
Interpolate to a resolution of 100m
calculate latitude and longitude,
convert to nc file

20211011
"""

import xarray as xr
import pandas as pd
import numpy as np
import subprocess
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import colors

def read_nc_point(filname,varname,ftime,lats,latn,lonl,lonr):
    f  = xr.open_dataset(filname) #yc,xc
    lat = f['latitude'].sel(xc=0).data
    lon = f['longitude'].sel(yc=0).data
    #ixc = f.xc[(lon>=lonl) & (lon<=lonr)]
    #iyc = f.yc[(lat>=lats) & (lat<=latn)]
    ixc = np.argmin(abs(lon-lonl))
    iyc = np.argmin(abs(lat-lats))
    print(ixc)
    print(iyc)
    if varname == "xwnd":
        u = f[varname].sel(time=ftime,yc=iyc,xc=ixc)
        v = f["ywnd"].sel(time=ftime,yc=iyc,xc=ixc)
        term = u
        term.data = np.power((np.power(u.data,2)+np.power(v.data,2)),0.5) # speed of wind
        var = term#.mean(dim=('yc','xc')) 
        del u, v, term
    else:
        #var = f[varname].sel(time=ftime,yc=iyc,xc=ixc).mean(dim=('yc','xc'))
        var = f[varname].sel(time=ftime,yc=iyc,xc=ixc)
    return var

# constants
BIGFONT=22
MIDFONT=18
SMFONT=14

lat = 20.95272
lon = 111.63450
lats = lat#-0.01 #21.7
latn = lat#+0.01 #21.9
lonl = lon#-0.02 #114.4
lonr = lon#+0.02 #114.6

filname = ["hsig","tps","wind"]
drawvar = ["hs"  ,"tps","xwnd"]
unit    = ["m"   ,"s"  ,"m/s" ]
nv = 0 

# define date of start time and forcast time
ftime = pd.date_range(start="2021-11-01 00",end="2021-11-10 00",freq="1H",closed=None)
case  = ["2021103012","2021103100", "2021103112"]
path  = "/home/metctm1/array74/Njord_Calypso/archive/calypso-gfs/"
figdir= "/home/lzhenn/cooperate/fig/%s.jpg"%(filname[nv])
title = "calypso-cfs %s (%s) (%.2fN,%.2fE) "%(filname[nv],unit[nv],lat,lon)

var = np.empty((len(case),len(ftime)), dtype = float)
for nc in range(len(case)):
    if nc >= 1: 
        filedir = path+case[nc]+"/calypso_"+filname[nv]+"_d01.nc"
    else:
        filedir = path+case[nc]+"/GBA_"+filname[nv]+"_d01.nc"
    var[nc,:] = read_nc_point(filedir,drawvar[nv],ftime,lats,latn,lonl,lonr) 
print(var)

fig = plt.figure(figsize=(12,9),dpi=300)
axe = fig.add_axes([0.05, 0.05, 0.9, 0.45])
axe.set_title(title,fontsize=MIDFONT)

axe.plot(ftime, np.mean(var,axis=0),color="blue")
std=np.std(var,axis=0)
axe.fill_between(ftime, np.mean(var,axis=0)-std,\
    np.mean(var,axis=0)+std,color="blue",alpha=0.25)

axe.set_ylabel("%s(%s)"%(filname[nv],unit[nv]),fontsize=MIDFONT)  # Add an x-label to the axes.
axe.set_xlabel("time",fontsize=MIDFONT)  # Add a y-label to the axes.
axe.tick_params(axis="both", which="major", labelsize=SMFONT)
#axe.xaxis.set_major_locator(mdates.MonthLocator(interval=4))
axe.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d_%H:00"))
plt.setp(axe.get_xticklabels(), rotation=30, ha="right")
#axe.legend(fontsize=MIDFONT)  # Add a legend.
fig.savefig(figdir,bbox_inches="tight",pad_inches=0.1)
#fig.savefig("%s%s.jpg"%(figdir,title),bbox_inches="tight",pad_inches=0.1)

