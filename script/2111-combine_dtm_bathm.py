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
filname=["/home/lzhenn/cooperate/data/Whole_HK_DTM_5m.nc",\
         "/home/lzhenn/cooperate/data/Whole_HK_DTM_100m.nc",\
         "/home/metctm1/array/data/Calypso/roms_d03.nc"]
fileout="/home/lzhenn/cooperate/data/dtm_bathymetry_100m3.nc"

f2  = xr.open_dataset(filname[2])
lats = f2['lat_rho'].sel(eta_rho=0 ,xi_rho=1).data
latn = f2['lat_rho'].sel(eta_rho=-1,xi_rho=1).data
lonl = f2['lon_rho'].sel(eta_rho=1,xi_rho=0).data
lonr = f2['lon_rho'].sel(eta_rho=1,xi_rho=-1).data
ilat = f2['lat_rho'].sel(xi_rho=1).data
ilon = f2['lon_rho'].sel(eta_rho=1).data

def main_run():
    #draw_jpg_map()
    #shengzhen_srtm('%s/map300dpi_100m_shengzhen.bmp'%figdir)
    combine_dtm_bath('%s/dtm_100m_runway.bmp'%figdir)
    ''' 
    ds = xr.open_dataset(fileout)
    var = ds['mask_rho'].data
    lat = ds['lat_rho'].data
    lon = ds['lon_rho'].data
    draw_map_bmp('%s/dtm_100m.jpg'%figdir,var,lat,lon)
    '''

def shengzhen_srtm(figname):
    var1 = f2['mask_rho']
    var = bmp2imask(figname,var1)
    del var1
    
    fstrm = xr.open_dataset('/home/lzhenn/cooperate/data/SRTM_N22E113114.nc')
    term = fstrm['hgt']
    print(term)
    hgt = term.interp(lat=ilat,lon=ilon,method='linear')
    print(hgt)
    hgt.data = np.ma.array(hgt.data, mask=var==0)
    print(hgt)
    ds = hgt.to_dataset(name="hgt")
    ds.to_netcdf('/home/lzhenn/cooperate/data/SRTM_shengzhen.nc')

def combine_dtm_bath(figname):
    # convert bmp2imask, then combine dtm and bathymetry
    # draw landsea_mask to test the combine results
    var1 = f2['mask_rho']

    f1  = xr.open_dataset(filname[1])
    dtm = f1['dtm'].sel(lat=ilat,lon=ilon,method="nearest")
    print(dtm)
    f3  = xr.open_dataset('/home/lzhenn/cooperate/data/SRTM_shengzhen.nc')
    term = f3['hgt'].data
    print(term)
    term = np.nan_to_num(term,nan=0)
    print(term)
    runway = bmp2imask('%s/only_runway2.bmp'%figdir,var1)
    runway = np.where(runway==1, 6, 0)
    print(runway)
    dtm.data = dtm.data + term + runway

    var = bmp2imask(figname,var1)
    var.attrs["option_0"] = "water"
    var.attrs["option_1"] = "land"
    ds = var.to_dataset(name="mask_rho")
    print(var)
    h = f2['h']
    h.values = np.where(var==0, -1*h, dtm)
    h.attrs["option<0"] = "water"
    h.attrs["option>0"] = "land"
    print(h)
    ds['h'] = h

    ds['lat_rho'] = f2['lat_rho']
    ds['lon_rho'] = f2['lon_rho']
    ds.to_netcdf(fileout)

def draw_jpg_map():
    # output the map to handle the data
    f1  = xr.open_dataset(filname[1])
    #f1  = xr.open_dataset(filname[0])
    var0 = f1['dtm'].data#.sel(lat=ilat,lon=ilon,method="nearest")
    ilat = f1['lat'].data
    ilon = f1['lon'].data
    var  = np.where(var0>0, 1, 0) # land 1, water 0
    var  = np.ma.array(var,mask=var==0)
    draw_map_bmp(figdir+"/map300dpi_5m.jpg",var,ilat,ilon)

def imask2bmp(var):
    im_w = var.shape[1]
    im_h = var.shape[0]
    
    image = Image.new('1', (im_w, im_h), 0)
    draw = ImageDraw.Draw(image)
    for x in range(im_w):
        for y in range(im_h):
            draw.point((x, y), fill=int(var[im_h-y-1,x]))
    image.save(figdir+'/bitmap_dtm.bmp', 'bmp')

def bmp2imask(figname,var):
    print("var.shape: ",var.shape)
    im_w = var.shape[1]
    im_h = var.shape[0]
    
    image = Image.open(figname)
    print(image.format, image.size, image.mode)
    image.thumbnail((im_w,im_h))
    values=list(image.getdata())
    for x in range(im_w):
        for y in range(im_h):
            var.values[im_h-y-1,x]=values[y*im_w+x]
    
    var.values = abs(var.values/255-1)
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
    #axe.set_xlim([ilon[0,0],ilon[-1,-1]])
    #axe.set_ylim([ilat[0,0],ilat[-1,-1]])
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
    #cont = axe.contourf(ilon, ilat, var, [0.5,2],zorder=0, 
    #         transform=ccrs.PlateCarree(),cmap=ncmap,extend="both",norm=norm1)

    axe.set_xlim([lonl,lonr])
    axe.set_ylim([lats,latn])
    #axe.set_xlim([ilon[0,0],ilon[-1,-1]])
    #axe.set_ylim([ilat[0,0],ilat[-1,-1]])
    axe.set_xticks(np.arange(np.ceil(ilon[0,0]),np.ceil(ilon[-1,-1]),0.5), crs=ccrs.PlateCarree())
    axe.set_yticks(np.arange(np.ceil(ilat[0,0]),np.ceil(ilat[-1,-1]),0.5), crs=ccrs.PlateCarree())
    axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
    axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))
    plt.savefig(figname,bbox_inches="tight",pad_inches=0.)
    del fig, axe

if __name__=='__main__':
    main_run()
