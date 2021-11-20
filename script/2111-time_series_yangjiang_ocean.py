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

def read_romsnc_point(filname,varname,ftime,lats,lonl):
    f  = xr.open_dataset(filname) #yc,xc
    if varname in ["speed","dir"]:
        lat = f['lat_u'].sel(xi_u=0).data
        lon = f['lon_u'].sel(eta_u=0).data
        ixc = np.argmin(abs(lon-lonl))
        iyc = np.argmin(abs(lat-lats))
        #u = f['u'].sel(ocean_time=ftime,s_rho=-0.016667,eta_u=iyc,xi_u=ixc)
        u = f['u'][0:len(ftime),29,iyc,ixc]
        
        lat = f['lat_v'].sel(xi_v=0).data
        lon = f['lon_v'].sel(eta_v=0).data
        ixc = np.argmin(abs(lon-lonl))
        iyc = np.argmin(abs(lat-lats))
        #v = f['v'].sel(ocean_time=ftime,s_rho=-0.016667,eta_v=iyc,xi_v=ixc)
        v = f['v'][0:len(ftime),29,iyc,ixc]
        
        if varname == "speed":
            var = np.power((np.power(u.data,2)+np.power(v.data,2)),0.5)*100 # speed of wind
        else:
            var = np.mod(180+np.rad2deg(np.arctan2(u,v)),360) # direction
            # 0 mean southerly wind
        del u, v
    return var

def read_nc_point(filname,varname,ftime,lats,lonl):
    f  = xr.open_dataset(filname) #yc,xc
    lat = f['latitude'].sel(xc=0).data
    lon = f['longitude'].sel(yc=0).data
    #ixc = f.xc[(lon>=lonl) & (lon<=lonr)]
    #iyc = f.yc[(lat>=lats) & (lat<=latn)]
    ixc = np.argmin(abs(lon-lonl))
    iyc = np.argmin(abs(lat-lats))
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
MIDFONT=20
SMFONT=18

lat = 20.961667 
lon = 111.608889
filname = ["hsig","tps","njord_his","njord_his"]
drawvar = ["hs"  ,"tps","speed","dir"]
unit    = ["m"   ,"s"  ,"cm/s" ,"deg"]
obsvar  = ["H3","Tp" ,'表  层','Unnamed: 14']

# define date of start time and forcast time
ftime = pd.date_range(start="2020-07-01 00",end="2020-07-11 00",freq="1H",closed=None)
path  = ["/home/lzhenn/cooperate/data/case_study/yangjiang-windfarm/2020063000_org/",
         "/home/lzhenn/array74/data/yangjiang-windfarm/report_expand.xlsx"]
case  = ["Coupled","Obs"]

for nv in range(0,4,1):
    figdir= "/home/lzhenn/cooperate/fig/%s.jpg"%(drawvar[nv])
    title = "ocean %s (%s) (%.2fN,%.2fE) "%(drawvar[nv],unit[nv],lat,lon)

    var = np.empty((len(case),len(ftime)), dtype = float)
    df = pd.read_excel(path[1], index_col=0, usecols=['日期',obsvar[nv]])
    var[1,:] = df[df.index.strftime('%Y-%m-%d %H').isin(ftime.strftime('%Y-%m-%d %H'))].to_numpy()[:,0]

    nc = 0
    for nt in range(0,len(ftime)-1,24):
        print(nt,ftime[nt])
        filedir = path[nc]+filname[nv]+"_d01."+ftime[nt].strftime("%Y%m%d")+".nc"
        if nv <= 1:
            var[nc,nt:(nt+25)] = read_nc_point(filedir,drawvar[nv],ftime[nt:(nt+25)],lat,lon) 
        else:
            var[nc,nt:(nt+25)] = read_romsnc_point(filedir,drawvar[nv],ftime[nt:(nt+25)],lat,lon) 

    fig = plt.figure(figsize=(12,9),dpi=300)
    axe = fig.add_axes([0.05, 0.05, 0.9, 0.45])
    axe.set_title(title,fontsize=MIDFONT) 
    for nc in range(len(var)):
        if nc == 1:
            axe.plot(ftime, var[nc,:], label=case[nc])
        else:
            rmse = np.sqrt(np.power(np.nan_to_num((var[nc,:] - var[1,:]),copy=False,nan=0),2).mean())
            mae = np.nan_to_num((var[nc,:] - var[1,:]),copy=False,nan=0).mean()
            axe.plot(ftime, var[nc,:], label="%s (RMSE: %.2f, MAE: %.2f)"%(case[nc],rmse,mae))
    #axe.set_ylabel("wind speed (m/s)",fontsize=MIDFONT)  # Add an x-label to the axes.
    axe.set_xlabel("time",fontsize=MIDFONT)  # Add a y-label to the axes.
    axe.tick_params(axis='both', which='major', labelsize=SMFONT)
    axe.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d_%H:00"))
    plt.setp(axe.get_xticklabels(), rotation=15, ha="right")
    axe.legend(fontsize=MIDFONT)  # Add a legend.
    fig.savefig(figdir,bbox_inches="tight",pad_inches=0.1)

