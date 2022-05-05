import xarray as xr
import pandas as pd
import numpy as np
import os, subprocess, wrf
from multiprocessing import Pool
from netCDF4 import Dataset
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.io.shapereader import Reader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmaps

title_font=18
label_font=14
plt.rcParams["font.weight"] = "bold"
font = {'family': 'sans-serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black',
        }

indir = ['/home/lzhenn/array74/data/ssp_rainfall/hist_2010s/',
         '/home/metctm1/array/data/imerg/daily/3B-DAY.MS.MRG.3IMERG.',
         '/home/metctm1/array/data/imerg/mon/3B-MO.MS.MRG.3IMERG.',
         '/home/metctm1/array/data/imerg/3hr']
case = ['hist-2010s','IMERG','IMERG','IMERG']
years = range(2011,2021)

def main_run():
    lat_sp = [2, 1, 0.5] 
    lon_sp = [2, 1, 0.5] 
    dom = ['d02','d03','d04']
    perid = ['Apr-Sep','AMJ','JAS']
    
    for nd in [0]:
        for ns in [0,1,2]:
            draw_precip_map(dom[nd],perid[ns],lat_sp[nd],lon_sp[nd])

def draw_precip_map(dom,season,lat_sp,lon_sp):
    if season == 'AMJ':
        datafile = ["%s/wrfout_%s.precc.*0[456].nc"%(indir[0],dom),
                    "%s20[1-2][0-9]0[456][0-3][0-9]-S000000-E235959.V06.nc4.nc4"%(indir[1]),
                    "%s20[1-2][1-9]0[456]01-S000000-E235959.0[456].V06B.HDF5.nc4"%(indir[2])]
        nday = 91
    if season == 'JAS':
        datafile = ["%s/wrfout_%s.precc.*0[789].nc"%(indir[0],dom),
                    "%s20[1-2][0-9]0[789][0-3][0-9]-S000000-E235959.V06.nc4.nc4"%(indir[1]),
                    "%s20[1-2][1-9]0[789]01-S000000-E235959.0[789].V06B.HDF5.nc4"%(indir[2])]
        nday = 92
    if season == 'Apr-Sep':
        datafile = ["%s/wrfout_%s.precc.*0[4-9].nc"%(indir[0],dom),
                    "%s20[1-2][0-9]0[4-9][0-3][0-9]-S000000-E235959.V06.nc4.nc4"%(indir[1]),
                    "%s20[1-2][1-9]0[4-9]01-S000000-E235959.0[4-9].V06B.HDF5.nc4"%(indir[2])]
        nday = 183
    
    bmlo = 0.45
    cnlevels = np.arange(10,1710,100) 
    fcolors = cmaps.precip2_17lev
    norm = colors.BoundaryNorm(boundaries=cnlevels, 
        ncolors=fcolors.N,extend='both')
    
    fig = plt.figure(figsize=(12,9),dpi=100)
    axe = fig.subplots(1,2, subplot_kw=dict(projection=ccrs.PlateCarree()))
    naxe = 0
    for nc in [0,2]:
        ax = axe[naxe]
        naxe = naxe+1
        ax.set_title("%s %s %s"%(dom,case[nc],season),fontsize=title_font,fontdict=font)
        
        coast_shp1 = Reader('/home/lzhenn/njord_pipeline/postprocess/shp/cnhimap.dbf').geometries()
        coastline1 = cfeat.ShapelyFeature(coast_shp1, ccrs.PlateCarree(), edgecolor="black", facecolor="none")
        ax.add_feature(coastline1, linewidth=0.8,zorder=1)
        if dom == 'd03':
            coast_shp2 = Reader('/home/lzhenn/njord_pipeline/postprocess/shp/gadm36_CHN_2.dbf').geometries()
            coastline2 = cfeat.ShapelyFeature(coast_shp2, ccrs.PlateCarree(), edgecolor="black", facecolor="none")
            ax.add_feature(coastline2, linewidth=0.8,zorder=1)
        
        fn_stream = subprocess.check_output('ls '+datafile[nc], shell=True).decode('utf-8')
        print(fn_stream)
        ds = xr.open_mfdataset(datafile[nc],concat_dim='time', combine='nested')
        if nc == 0:
            var = ds['temp']
            lats = ds['lat'][0]
            latn = ds['lat'][-1]
            lonl = ds['lon'][0]
            lonr = ds['lon'][-1]
        else:
            lat = ds.lat
            lon = ds.lon
            ilon = lon[(lon>=lonl) & (lon<=lonr)]
            ilat = lat[(lat>=lats) & (lat<=latn)]
            var = ds['precipitation'].sel(lat=ilat,lon=ilon)
            #var = ds['HQprecipitation'].sel(lat=ilat,lon=ilon)
            var = var.transpose('time','lat','lon')
        var = var.mean('time')
        print(case[nc])
        print(var)
        cont = ax.contourf(var.lon, var.lat, var*24*nday, cnlevels, 
            transform=ccrs.PlateCarree(),cmap=fcolors, extend='both',norm=norm)
        
        ax.set_xticks(np.arange(np.ceil(lonl),np.ceil(lonr),lon_sp), crs=ccrs.PlateCarree())
        ax.set_yticks(np.arange(np.ceil(lats),np.ceil(latn),lat_sp), crs=ccrs.PlateCarree())
        ax.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
        ax.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))
        
    position = fig.add_axes([0.2, bmlo+0.005, 0.7, 0.01]) #left, bottom, width, height
    cb = plt.colorbar(cont, cax=position ,orientation='horizontal')#, shrink=.9)
    cb.set_label(label='precip (mm/month)', size=title_font) #, weight='bold'

    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    fig.savefig("/home/lzhenn/cooperate/fig/obs_cmip6_%s_%s"%(dom,season)
        ,bbox_inches='tight',pad_inches=0.1)

if __name__=='__main__':
    main_run()
