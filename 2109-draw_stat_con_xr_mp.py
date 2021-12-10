#!/usr/bin/env python
import sys
import subprocess
import xarray as xr
import numpy as np
import gc #garbage collector
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmaps
#matplotlib.use('Agg')

lonl=0  #0  #
lonr=150#360#
lats=15 #10 #
latn=70 #90 #
lat_sp = 20 #30
lon_sp = 30 #60
nrow = 3
ncol = 3
bmlo = 0.45 #0.55 #
title_font=14
label_font=10
xbar=[0.05,0.37,0.69]

# level == 0 #small range
cnlvl=[[0    ,320  ,20  ], # 0Feature Density
       [0    ,3.2  ,0.2 ], # 1Genesis Density
       [0    ,3.2  ,0.2 ], # 2Lysis Density
       [-1   ,1    ,1   ], # 3Mean Area
       [-0.8 ,0.8  ,0.1 ], # 4Mean Growth/Decay Rate
       [-1   ,1    ,1   ], # 5Mean Anisotropy
       [0    ,8    ,0.5 ], # 6Mean Lifetime
       [0    ,80   ,5   ], # 7Mean Speed
       [0    ,16   ,1   ], # 8Mean Intensity
       [-1.6 ,1.6  ,0.2 ], # 9Mean Tendency
       [-1   ,1    ,1   ], # 10Spare1
       [-1   ,1    ,1   ], # 11Spare2
       [0    ,80   ,5   ], # 12Std of Speed
       [0    ,3.2  ,0.2 ], # 13Std of Intensity
       [0    ,16   ,1   ], # 14Track Density
       [-1   ,1    ,1   ], # 15X-component of Mean Orientation Vector
       [-40  ,40   ,5   ], # 16X-component of Mean Velocity
       [-1   ,1    ,1   ], # 17Y-component of Mean Orientation Vector
       [-40  ,40   ,5  ]] # 18Y-component of Mean Velocity

draw_var = ["fden","gden","lden","marea","mgdr","",
            "mlif","msp" ,"mstr","mten" ,""    ,"",
            ""    ,""    ,"tden",""    ] # 7 variables
#draw=[8,9,6]
draw=[1,2,14]
#draw=[14]
#draw=[1,2,14,8,9,6]
lev = [250,500,850]
if len(sys.argv) < 2 :
    prefix = "ff"
    suffix = ''
    dbox = 0 
    level = 2
else:
    prefix = sys.argv[1]  #'ff_250_500_no'
    suffix = sys.argv[2]  #'ff_250_500_no'
    level = int(sys.argv[3])
    dbox  = int(sys.argv[4])

if dbox >= 1 :
    flats = int(sys.argv[5])
    flatn = int(sys.argv[6])
    flonl = int(sys.argv[7])
    flonr = int(sys.argv[8])

if level == 1: #middle range 
    cnlvl[ 1][:]=[0    ,4.8  ,0.3 ]
    cnlvl[ 2][:]=[0    ,4.8  ,0.3 ]
    cnlvl[14][:]=[0    ,24   ,1.5 ]
if level == 2: # large range,use total cnlvlel bar 
    cnlvl[ 0][:]=[0    ,1600 ,100 ]
    cnlvl[ 1][:]=[0    ,8    ,0.5 ]
    cnlvl[ 2][:]=[0    ,8    ,0.5 ]
    cnlvl[14][:]=[0    ,32   ,2   ]

files = '/home/users/qd201969/ERA5-1HR-lev/statistic/ff_250_1980-2020_stat.nc'
f = xr.open_dataset(files)
lat = f.lat
lon = f.long
ilon = lon[(lon>=lonl) & (lon<=lonr)]
ilat = lat[(lat>=lats) & (lat<=latn)]
del f,lat,lon

ds = xr.open_dataset('/home/users/qd201969/data/ERA5_mon_u_1979-2020.nc')
da = ds['u'].sel(level=200,longitude=ilon,latitude=ilat,method="nearest").load()
# increased performance by loading data into memory first, e.g., with load()
uwnd = da.mean('time')
del ds, da

ds = xr.open_dataset("/home/users/qd201969/gtopo30_0.9x1.25.nc")
phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").load()
phis = phis/9.8 # transfer from m2/s2 to m
del ds
gc.collect()

figdir = "/home/users/qd201969/uor_track/fig/stat_annual"+suffix
fig = plt.figure(figsize=(12,12),dpi=300)
ax = fig.subplots(nrow, ncol, subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=180.0))) #sharex=True, sharey=True
for nl in range(0,len(lev),1):
    files = '/home/users/qd201969/ERA5-1HR-lev/statistic/%s_%d_1980-2020%s_stat.nc'%(prefix,lev[nl],suffix)
    print(files)
    f = xr.open_dataset(files)
    for nv in range(0,len(draw),1):#,len(f),1):
        var1 = f[draw_var[draw[nv]]][0,1,2]
        var = f[draw_var[draw[nv]]].sel(long=ilon,lat=ilat).mean("time")
        print(var)
        if draw[nv] == 9:
            var=var*24
        if draw[nv] > 2 and draw[nv] != 14:
            tden = f['tden'].sel(long=ilon,lat=ilat).load()
            mask = tden < 1.0
            var.values=np.ma.array(var.values,mask=mask)
        
        cnlevels = np.arange(cnlvl[draw[nv]][0], cnlvl[draw[nv]][1]+cnlvl[draw[nv]][2], cnlvl[draw[nv]][2])
        if cnlvl[draw[nv]][0] < 0 :
            fcolors = cmaps.BlueDarkRed18
        else:
            fcolors = cmaps.precip2_17lev
        norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=fcolors.N,extend='both')
        
        axe = ax[nl][nv]
        axe.add_feature(cfeat.GSHHSFeature(levels=[1,2],edgecolor='k'), linewidth=0.8, zorder=1)
        axe.set_title("%d %s"%(lev[nl],var1.long_name),fontsize=title_font)

        shad = axe.contourf(ilon, ilat, var, cnlevels, 
                     transform=ccrs.PlateCarree(),cmap=fcolors,extend='both',norm=norm)
        topo = axe.contour(ilon, ilat, phis, [1500,3000,4500], 
                     transform=ccrs.PlateCarree(),colors='black',linewidths=1.5)
        if dbox >= 1 :
            axe.plot([flonl,flonl,flonr,flonr,flonl],[flatn,flats,flats,flatn,flatn], 
                     linewidth=2.5, color='black', transform=ccrs.PlateCarree()) # filter box

        jets = axe.contour(ilon, ilat, uwnd, [25,100], 
                     transform=ccrs.PlateCarree(),colors='darkviolet',linewidths=2)

        if nv == 0:
            axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
            axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
        if nl == (nrow-1):
            axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
            axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))
            position = fig.add_axes([xbar[nv], bmlo+0.05, 0.28, 0.008]) #left, bottom, width, height
            cb = plt.colorbar(shad, cax=position ,orientation='horizontal')#, shrink=.9)

plt.tight_layout(rect=(0,bmlo,1,1))
plt.savefig(figdir+".png", bbox_inches='tight',pad_inches=0.01)


