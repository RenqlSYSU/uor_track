import xarray as xr
import numpy as np
import subprocess
import os 
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.colors import ListedColormap
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.io.shapereader import Reader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.interpolate import griddata
from PIL import Image, ImageDraw 
import cmaps

title_font=18
label_font=14
plt.rcParams["font.weight"] = "bold"
font = {'family': 'sans-serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black',
        }

figdir = "/home/lzhenn/cooperate/fig"
path = '/home/lzhenn/cooperate/data/SRTM'
SAMPLES = 1201  # Change this to 3601 for SRTM1

def main_run():
    #save_hgt_nc()
    #fstrm = xr.open_dataset('/home/lzhenn/cooperate/data/SRTM_shengzhen.nc')
    #fstrm = xr.open_dataset('/home/lzhenn/cooperate/data/SRTM_N22E113114.nc')
    #hgt = fstrm['hgt']
    #draw_map('%s/SRTM_hgt.jpg'%figdir,'SRTM 1201',hgt,hgt.lat,hgt.lon)
    fileout="/home/lzhenn/cooperate/data/dtm_bathymetry_100m3.nc"
    f1  = xr.open_dataset(fileout)
    var = f1['h'].data
    ilat = f1['lat_rho'].sel(xi_rho=1).data
    ilon = f1['lon_rho'].sel(eta_rho=1).data
    draw_map('%s/bathymetry_roms_d03.jpg'%figdir,'h',var,ilat,ilon)

def save_hgt_nc():
    lats = 22
    latn = 22
    lonl = 113
    lonr = 114
    lat = np.arange( lats, latn+1/(SAMPLES-1)+1, 1/(SAMPLES-1) )
    lon = np.arange( lonl, lonr+1/(SAMPLES-1)+1, 1/(SAMPLES-1) )
    ele1 = read_hgt(lats,lonl)
    ele2 = read_hgt(lats,lonr)
    ele = np.concatenate((ele1,ele2[:,1:]),axis=1)
    draw_map('%s/SRTM_hgt.jpg'%figdir,'SRTM 1201',ele[::-1,:],lat,lon)

    ds = xr.Dataset(
             {
                     "hgt" : (["lat","lon"], ele[::-1,:]),
                     },
             coords={
                     "lat" : (["lat"], lat),
                     "lon" : (["lon"], lon),
                     },
             )
    ds.attrs["description"] = "Space resolution 90m. "+\
            "The data is computed from https://cgiarcsi.community/data/srtm-90m-digital-elevation-database-v4-1/?amp" 
    ds.to_netcdf("/home/lzhenn/cooperate/data/SRTM_N22E113114.nc","w")

    '''
    for i in range(lats,latn+1):
        for j in range(lonl,lonr+1):
            if os.path.exists(f_name):
                lat, lon, ele = read_hgt(f_name, i, j)
                draw_map('%s/SRTM_hgt.jpg'%figdir,'SRTM 1201',ele[::-1,:],lat,lon)
            else:
                print('%s do not exist'%f_name)
    '''
def read_hgt(lat, lon):
    f_name = '%s/N%dE%d.hgt'%(path,lat,lon)
    with open(f_name, 'rb') as hgt_data:
        elevations = np.fromfile(hgt_data, np.dtype('>i2'), SAMPLES * SAMPLES) \
            .reshape((SAMPLES, SAMPLES))
    #lat_range = np.arange(lat, lat + 1 / (SAMPLES-1) + 1, 1 / (SAMPLES-1))
    #lon_range = np.arange(lon, lon + 1 + 1 / (SAMPLES-1), 1 / (SAMPLES-1))
    return elevations

def draw_map(figname,title,var,ilat,ilon):
    fig = plt.figure(figsize=(12,9),dpi=300)
    axe = plt.axes(projection=ccrs.PlateCarree())
    
    axe.set_title(title, fontsize=title_font, fontdict=font)
    coast_shp2 = Reader(os.getenv("SHP_LIB")+"/china_coast/china_coastline.dbf").geometries()
    coastline2 = cfeat.ShapelyFeature(coast_shp2, ccrs.PlateCarree(), edgecolor="red",facecolor="none")
    axe.add_feature(coastline2, linewidth=0.8,zorder=2)
    
    cnlevels = np.arange(0,17,1) 
    ncmap = cmaps.precip2_17lev
    norm  = colors.BoundaryNorm(boundaries=cnlevels, ncolors=ncmap.N,extend="both")
    cont = axe.pcolormesh(ilon, ilat, var, zorder=0, 
             transform=ccrs.PlateCarree(),cmap=ncmap,norm=norm)
    plt.colorbar(cont, ax=axe, shrink=.7)

    axe.set_xlim([ilon[0],ilon[-1]])
    axe.set_ylim([ilat[0],ilat[-1]])
    axe.set_xticks(np.arange(np.ceil(ilon[0]),np.ceil(ilon[-1]),0.5), crs=ccrs.PlateCarree())
    axe.set_yticks(np.arange(np.ceil(ilat[0]),np.ceil(ilat[-1]),0.5), crs=ccrs.PlateCarree())
    axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
    axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))
    plt.savefig(figname,bbox_inches="tight",pad_inches=0.)

if __name__=='__main__':
    main_run()
