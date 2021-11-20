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

def imask2bmp(var):
    im_w = var.shape[1]
    im_h = var.shape[0]
    
    image = Image.new('1', (im_w, im_h), 0)
    draw = ImageDraw.Draw(image)
    for x in range(im_w):
        for y in range(im_h):
            draw.point((x, y), fill=int(var[im_h-y-1,x]))
    image.save(figdir+'bitmap_dtm.bmp', 'bmp')

def fig_resize(figname,var):
    print("var.shape: ",var.shape)
    im_w = var.shape[1]
    im_h = var.shape[0]
    image = Image.open(figname)
    print(image.format, image.size, image.mode)
    image.thumbnail((im_w,im_h))
    image.save(figname,'bmp')
   
def bmp2imask(figname,var):
    print("var.shape: ",var.shape)
    im_w = var.shape[1]
    im_h = var.shape[0]
    
    image = Image.open(figname)
    print(image.format, image.size, image.mode)
    values=list(image.getdata())
    for x in range(im_w):
        for y in range(im_h):
            var.values[im_h-y-1,x]=values[y*im_w+x]
    
    var.values = abs(var.values/255-1)
    return var

def draw_map(figname,var,fmap=False):
    fig = plt.figure(figsize=(12,9),dpi=300)
    axe = plt.axes(projection=ccrs.PlateCarree())
   
    if fmap == True:
        coast_shp  = Reader(os.getenv("SHP_LIB")+"/china_coast/china_coastline.dbf").geometries()
        coastline1 = cfeat.ShapelyFeature(coast_shp, ccrs.PlateCarree(), edgecolor="black", facecolor="black")
        axe.add_feature(coastline1, linewidth=0.8,zorder=1)
    
    #coast_shp2 = Reader(os.getenv("SHP_LIB")+"/china_coast/china_coastline.dbf").geometries()
    #coastline2 = cfeat.ShapelyFeature(coast_shp2, ccrs.PlateCarree(), edgecolor="red",facecolor="none")
    #axe.add_feature(coastline2, linewidth=0.8,zorder=2)
    
    ncmap = ListedColormap(["white","black","green"])
    norm1 = colors.BoundaryNorm(boundaries=[0.5,2], ncolors=3,extend="both")
    cont = axe.pcolormesh(ilon, ilat, var, zorder=0, 
             transform=ccrs.PlateCarree(),cmap=ncmap,norm=norm1)
    #cont = axe.contourf(ilon, ilat, var, [0.5,2],zorder=0, 
    #         transform=ccrs.PlateCarree(),cmap=ncmap,extend="both",norm=norm1)

    axe.set_xlim([lonl,lonr])
    axe.set_ylim([lats,latn])
    plt.axis('off')
    plt.savefig(figname,bbox_inches="tight",pad_inches=0.)
    del fig, axe

figdir = "/home/lzhenn/cooperate/fig/"
filname=["/home/lzhenn/cooperate/data/Whole_HK_DTM_5m.nc",\
         "/home/lzhenn/cooperate/data/Whole_HK_DTM_100m.nc",\
         "/home/metctm1/array/data/Calypso/roms_d03.nc"]
fileout="/home/lzhenn/cooperate/data/dtm_bathymetry_100m.nc"

f2  = xr.open_dataset(filname[2])
lats = f2['lat_rho'].sel(eta_rho=0 ,xi_rho=1).data
latn = f2['lat_rho'].sel(eta_rho=-1,xi_rho=1).data
lonl = f2['lon_rho'].sel(eta_rho=1,xi_rho=0).data
lonr = f2['lon_rho'].sel(eta_rho=1,xi_rho=-1).data

ilat = f2['lat_rho'].sel(xi_rho=1).data
ilon = f2['lon_rho'].sel(eta_rho=1).data

step = 2
if step == 0:
    # output the map to handle the data
    #f1  = xr.open_dataset(filname[1])
    f1  = xr.open_dataset(filname[0])
    var0 = f1['dtm']#.sel(lat=ilat,lon=ilon,method="nearest")
    ilat = f1['lat']
    ilon = f1['lon']
    var  = np.where(var0>0, 1, 0) # land 1, water 0
    draw_map(figdir+"map300dpi_5m.jpg",var,fmap=True)

figname = figdir+"map300dpi_5m_remove.bmp" 
if step == 1:
    var1 = f2['mask_rho']
    var = fig_resize(figname,var1)

step = 2
if step == 2:
    # convert bmp2imask, then combine dtm and bathymetry
    # draw landsea_mask to test the combine results
    var1 = f2['mask_rho']
    var = bmp2imask(figname,var1)
    del var1
    var.attrs["option_0"] = "water"
    var.attrs["option_1"] = "land"
    ds = var.to_dataset(name="mask_rho")
    print(var)
    draw_map(figdir+"test_landsea_mask.jpg",var)

    f1  = xr.open_dataset(filname[1])
    var0 = f1['dtm'].sel(lat=ilat,lon=ilon,method="nearest")
    h = f2['h']

    h.values = np.where(var==0, -1*h, var0)
    h.attrs["option<0"] = "water"
    h.attrs["option>0"] = "land"
    print(h)
    ds['h'] = h

    ds['lat_rho'] = f2['lat_rho']
    ds['lon_rho'] = f2['lon_rho']
    ds.to_netcdf(fileout)

