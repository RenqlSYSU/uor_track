"""
read Whole_HK_DTM_6m.asc, 
Interpolate to a resolution of 100m
calculate latitude and longitude,
convert to nc file

20211011
"""

import xarray as xr
import pandas as pd
import numpy as np
from datetime import datetime
import subprocess
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import colors
import cmaps

def read_romsnc_point(filname,varname,ftime,lats,lonl):
    ds = xr.open_dataset(filname) #yc,xc
    if varname in ["temp"]:
        lat = ds['lat_rho'].sel(xi_rho=0).data
        lon = ds['lon_rho'].sel(eta_rho=0).data
        ixc = np.argmin(abs(lon-lonl))
        iyc = np.argmin(abs(lat-lats))
        var = ds[varname].sel(ocean_time=ftime,eta_rho=iyc,xi_rho=ixc)
        
        zeta = ds.zeta.sel(ocean_time=ftime,eta_rho=iyc,xi_rho=ixc)
        h = ds.h.sel(eta_rho=iyc,xi_rho=ixc)
        if ds.Vtransform == 1:
            Zo_rho = ds.hc * (ds.s_rho - ds.Cs_r) + ds.Cs_r * h 
            z_rho = Zo_rho + zeta * (1 + Zo_rho / h)
        elif ds.Vtransform == 2:
            Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * h) / (ds.hc + h)
            z_rho = zeta + (zeta + h) * Zo_rho
        print(z_rho[0,:])

        depth = np.arange(-100,2,2)
        data = np.empty((len(ftime),len(depth)), dtype = float)
        for nt in range(len(ftime)):
            var.coords['s_rho'] = z_rho[nt,:].data
            data[nt,:] = var.interp(ocean_time=ftime[nt], s_rho=depth, method="cubic").data
        var = xr.DataArray(data, coords=[("ocean_time", ftime), ("z_rho", depth)])
        
    if varname in ["speed","dir"]:
        lat = ds['lat_u'].sel(xi_u=0).data
        lon = ds['lon_u'].sel(eta_u=0).data
        ixc = np.argmin(abs(lon-lonl))
        iyc = np.argmin(abs(lat-lats))
        #u = ds['u'].sel(ocean_time=ftime,s_rho=-0.016667,eta_u=iyc,xi_u=ixc)
        u = ds['u'][0:len(ftime),29,iyc,ixc]
        
        lat = ds['lat_v'].sel(xi_v=0).data
        lon = ds['lon_v'].sel(eta_v=0).data
        ixc = np.argmin(abs(lon-lonl))
        iyc = np.argmin(abs(lat-lats))
        #v = ds['v'].sel(ocean_time=ftime,s_rho=-0.016667,eta_v=iyc,xi_v=ixc)
        v = ds['v'][0:len(ftime),29,iyc,ixc]
        
        if varname == "speed":
            var = np.power((np.power(u.data,2)+np.power(v.data,2)),0.5)*100 # speed of wind
        else:
            var = np.mod(180+np.rad2deg(np.arctan2(u,v)),360) # direction
            # 0 mean southerly wind
        del u, v
    return var

# constants
MIDFONT=18
SMFONT=12

lat = 19.236
lon = 118.411
drawvar = ["temp"]
unit    = ["Celsius"]
nv = 0

# define date of start time and forcast time
ftime = pd.date_range(start="2018-09-14 12",end="2018-09-16 12",freq="0.5H",closed=None)
path  = "/home/metctm1/array/data/1911-COAWST/ERA5_TY2001_org/gba_ocean_his.nc" 
figdir= "/home/lzhenn/cooperate/fig/roms_%s.jpg"%(drawvar[nv])
title = "ocean %s (%s) (%.2fN,%.2fE) "%(drawvar[nv],unit[nv],lat,lon)

var = read_romsnc_point(path, drawvar[nv], ftime, lat, lon)
print(var)

fig = plt.figure(figsize=(12,9),dpi=300)
axe = fig.add_axes([0.05, 0.05, 0.9, 0.4])
axe.set_title(title,fontsize=MIDFONT) 

cnlevels = np.arange(20,29.5,0.25) #500Z, gpm
ncmap = colors.ListedColormap(cmaps.MPL_jet(range(0,127,3))) #cmaps.precip2_17lev
norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=ncmap.N,extend='both')

cont = axe.contourf(var.ocean_time, var.z_rho, var.transpose(), cnlevels, 
     cmap=ncmap,extend='both',norm=norm)
plt.colorbar(cont, ax=axe)

refer = datetime(2018, 9, 15, 12, 00)
axe.plot([refer,refer],[-100,0],linewidth=2.0,color='black',linestyle='dashed')

axe.set_ylim([-100,0])
axe.set_ylabel("depth (m)",fontsize=MIDFONT)  # Add an x-label to the axes.
#axe.set_xlabel("time",fontsize=MIDFONT)  # Add a y-label to the axes.
axe.tick_params(axis='both', which='major', labelsize=SMFONT)
axe.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d_%H:%M"))
plt.setp(axe.get_xticklabels(), rotation=15, ha="right")
fig.savefig(figdir,bbox_inches="tight",pad_inches=0.1)

