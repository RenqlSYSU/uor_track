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
from PIL import Image, ImageDraw 

figdir = "/home/lzhenn/cooperate/fig"
filname="/home/lzhenn/cooperate/data/dtm_bathymetry_100m3.nc"
romdata='/home/metctm1/array/data/Calypso/roms_d03.nc'
fileout="/home/lzhenn/cooperate/data/dtm_bathymetry_100m_lantau.nc"

f2  = xr.open_dataset(filname)
lats = f2['lat_rho'].sel(eta_rho=0 ,xi_rho=1).data
latn = f2['lat_rho'].sel(eta_rho=-1,xi_rho=1).data
lonl = f2['lon_rho'].sel(eta_rho=1,xi_rho=0).data
lonr = f2['lon_rho'].sel(eta_rho=1,xi_rho=-1).data
lat = f2['lat_rho'].data
lon = f2['lon_rho'].data

def main_run():
    #var = f2['mask_rho'].data
    #imask2bmp(var)
    #draw_map_bmp('%s/dtm.jpg'%figdir,var,lat,lon)
    
    combine_dtm_bath('%s/bitmap_dtm_lantau.bmp'%figdir)

def combine_dtm_bath(figname):
    # convert bmp2imask, obtain new landmask
    # change lantau height to 4m, obtain new height 
    # draw landsea_mask to test the combine results
    var1 = f2['mask_rho'].load()
    var = bmp2imask(figname,var1)
    print(var)
    ds = var.to_dataset(name="mask_rho")
    
    dtm = f2['h'].load()
    runway = bmp2imask('%s/only_lantau.bmp'%figdir,var1)
    dtm.data = np.where(runway.data==1, 4, dtm.data)
    
    f1 = xr.open_dataset(romdata)
    h = f1['h']
    dtm.data = np.where(var.data==0, -1*h, dtm.data)
    ds['h'] = dtm
    ds['lat_rho'] = f2['lat_rho']
    ds['lon_rho'] = f2['lon_rho']
    ds.to_netcdf(fileout)

    ds = xr.open_dataset(fileout)
    var = ds['mask_rho'].data
    draw_map_bmp('%s/dtm_100m.jpg'%figdir,var,lat,lon)

def imask2bmp(var):
    im_w = var.shape[1]
    im_h = var.shape[0]
    
    image = Image.new('1', (im_w, im_h), 0)
    draw = ImageDraw.Draw(image)
    for x in range(im_w):
        for y in range(im_h):
            draw.point((x, y), fill=int(abs(var[im_h-y-1,x]-1)))
    image.save(figdir+'/bitmap_dtm.bmp', 'bmp')

def bmp2imask(figname,term):
    var = term.copy()
    print("var.shape: ",var.shape)
    im_w = var.shape[1]
    im_h = var.shape[0]
    
    image = Image.open(figname)
    print(image.format, image.size, image.mode)
    image.thumbnail((im_w,im_h))
    values=list(image.getdata())
    for x in range(im_w):
        for y in range(im_h):
            var.data[im_h-y-1,x]=values[y*im_w+x]
    
    var.data = abs(var.data/255-1)
    return var

def draw_map_bmp(figname,var,ilat,ilon):
    fig = plt.figure(figsize=(12,9),dpi=300)
    axe = plt.axes(projection=ccrs.PlateCarree())
    coast_shp2 = Reader(os.getenv("SHP_LIB")+"/china_coast/china_coastline.dbf").geometries()
    coastline2 = cfeat.ShapelyFeature(coast_shp2, ccrs.PlateCarree(), edgecolor="red",facecolor="none")
    axe.add_feature(coastline2, linewidth=0.5,zorder=2)
    
    ncmap = ListedColormap(["white","black","green"])
    norm1 = colors.BoundaryNorm(boundaries=[0.5,2], ncolors=3,extend="both")
    cont = axe.pcolormesh(ilon, ilat, var, zorder=1, 
             transform=ccrs.PlateCarree(),cmap=ncmap,norm=norm1)

    axe.set_xlim([lonl,lonr])
    axe.set_ylim([lats,latn])
    plt.axis('off')
    plt.savefig(figname,bbox_inches="tight",pad_inches=0.)
    del fig, axe

def draw_map(figname,title,var,ilat,ilon):
    fig = plt.figure(figsize=(12,9),dpi=300)
    axe = plt.axes(projection=ccrs.PlateCarree())
    axe.set_title(title, fontsize=title_font, fontdict=font)
    coast_shp2 = Reader(os.getenv("SHP_LIB")+"/china_coast/china_coastline.dbf").geometries()
    coastline2 = cfeat.ShapelyFeature(coast_shp2, ccrs.PlateCarree(), edgecolor="red",facecolor="none")
    axe.add_feature(coastline2, linewidth=0.8,zorder=2)
    
    ncmap = ListedColormap(["white","black","green"])
    norm1 = colors.BoundaryNorm(boundaries=[0.5,2], ncolors=3,extend="both")
    cont = axe.pcolormesh(ilon, ilat, var, zorder=0, 
             transform=ccrs.PlateCarree(),cmap=ncmap,norm=norm1)

    axe.set_xlim([lonl,lonr])
    axe.set_ylim([lats,latn])
    axe.set_xticks(np.arange(np.ceil(ilon[0,0]),np.ceil(ilon[-1,-1]),0.5), crs=ccrs.PlateCarree())
    axe.set_yticks(np.arange(np.ceil(ilat[0,0]),np.ceil(ilat[-1,-1]),0.5), crs=ccrs.PlateCarree())
    axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
    axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))
    plt.savefig(figname,bbox_inches="tight",pad_inches=0.)
    del fig, axe

if __name__=='__main__':
    main_run()
