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
from scipy.io import loadmat
from scipy.interpolate import griddata
from PIL import Image, ImageDraw 

title_font=18
label_font=14
plt.rcParams["font.weight"] = "bold"
font = {'family': 'sans-serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black',
        }

figdir = "/home/lzhenn/cooperate/fig"
path = '/home/metctm1/array/data/Calypso'
def main_run():
    #land2romsd02_jpg()
    bmp2romsd02()

def bmp2romsd02():
    ds = xr.open_dataset('%s/roms_d02.nc'%path)
    dst_lat = ds['lat_rho'].data
    dst_lon = ds['lon_rho'].data
    var = ds['mask_rho'].data
    var_rho = bmp2imask('%s/interp_land300dpi_one.bmp'%figdir,var)
    draw_map_bmp('%s/roms_d02_land.png'%figdir,var_rho,dst_lat,dst_lon)
    #var_rho = np.where(var_rho==0, 1, 0)
    ds['mask_rho'].data = var_rho
    
    dst_lat2 = ds['lat_u'].data
    dst_lon2 = ds['lon_u'].data
    var = ds['mask_u'].data
    draw_map('%s/roms_d02_mask_u.png'%figdir,'roms_d02',var,dst_lat2,dst_lon2)
    dst_var2 = interp_griddata(var_rho,dst_lat,dst_lon,dst_lat2,dst_lon2)
    draw_map('%s/roms_d02_mask_u2.png'%figdir,'New roms_d02',dst_var2,dst_lat2,dst_lon2)
    ds['mask_u'].data = dst_var2
    print(dst_var2)
    
    dst_lat2 = ds['lat_v'].data
    dst_lon2 = ds['lon_v'].data
    var = ds['mask_v'].data
    draw_map('%s/roms_d02_mask_v.png'%figdir,'roms_d02',var,dst_lat2,dst_lon2)
    dst_var2 = interp_griddata(var_rho,dst_lat,dst_lon,dst_lat2,dst_lon2)
    draw_map('%s/roms_d02_mask_v2.png'%figdir,'New roms_d02',dst_var2,dst_lat2,dst_lon2)
    ds['mask_v'].data = dst_var2
    print(dst_var2)
    
    dst_lat2 = ds['lat_psi'].data
    dst_lon2 = ds['lon_psi'].data
    var = ds['mask_psi'].data
    draw_map('%s/roms_d02_mask_psi.png'%figdir,'roms_d02',var,dst_lat2,dst_lon2)
    dst_var2 = interp_griddata(var_rho,dst_lat,dst_lon,dst_lat2,dst_lon2)
    draw_map('%s/roms_d02_mask_psi2.png'%figdir,'New roms_d02',dst_var2,dst_lat2,dst_lon2)
    ds['mask_psi'].data = dst_var2
    print(dst_var2)
    
    ds.to_netcdf('/home/lzhenn/cooperate/data/new_roms_d02.nc','w')

def land2romsd02_jpg():
    land = loadmat('%s/land.mat'%path)['land']
    lat = loadmat('%s/lat.mat'%path)['lat']
    lon = loadmat('%s/lon.mat'%path)['lon']
    #imask2bmp('%s/land.bmp'%figdir,land) # 0 black
    #draw_map_bmp('%s/land600dpi.jpg'%figdir,land,lat,lon) # 1 black
    #draw_map('%s/land.png'%figdir,'High land mask',land,lat,lon)
    
    ds = xr.open_dataset('%s/roms_d02.nc'%path)
    dst_lat = ds['lat_rho'].data
    dst_lon = ds['lon_rho'].data
    dst_var = interp_griddata(land,lat,lon,dst_lat,dst_lon)
    draw_map_bmp('%s/interp_land300dpi.jpg'%figdir,dst_var,dst_lat,dst_lon) # 1 black
    #draw_map('%s/roms_d02_land2.png'%figdir,'New roms_d02',dst_var,dst_lat,dst_lon)
    print(dst_var)
    
def imask2bmp(figname,var):
    im_w = var.shape[1]
    im_h = var.shape[0]
    image = Image.new('1', (im_w, im_h), 0)
    draw = ImageDraw.Draw(image)
    for x in range(im_w):
        for y in range(im_h):
            draw.point((x, y), fill=int(var[im_h-y-1,x]))
    image.save(figname, 'bmp')

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
            var[im_h-y-1,x]=values[y*im_w+x]
    
    var = var/255
    print(var)
    print(var.dtype)
    return var

def interp_griddata(var,lat,lon,dst_lat,dst_lon):
    '''
    All input variables are 2D numpy array
    1. compress all the variables into one dimension
    2. use scipy.interpolate.griddata to interp 
    '''
    print('initial grid: ',lat.shape)
    print('dst_grid: ',dst_lat.shape)
    points = np.column_stack((lat.flatten(),lon.flatten())) #(nlat*nlon,2)
    dst_points = np.column_stack((dst_lat.flatten(),
                                  dst_lon.flatten())) #(nlat*nlon,2)
    dst_var = griddata(points, var.flatten(), dst_points,
        method='nearest')
    return dst_var.reshape(dst_lat.shape)

def draw_map_bmp(figname,var,ilat,ilon):
    fig = plt.figure(figsize=(12,9),dpi=300)
    axe = plt.axes(projection=ccrs.PlateCarree())
    coast_shp2 = Reader(os.getenv("SHP_LIB")+"/china_coast/china_coastline.dbf").geometries()
    coastline2 = cfeat.ShapelyFeature(coast_shp2, ccrs.PlateCarree(), edgecolor="red",facecolor="none")
    axe.add_feature(coastline2, linewidth=0.5,zorder=2)
    
    ncmap = ListedColormap(["white","black","green"])
    norm1 = colors.BoundaryNorm(boundaries=[0.5,2], ncolors=3,extend="both")
    cont = axe.pcolormesh(ilon, ilat, var, zorder=0, 
             transform=ccrs.PlateCarree(),cmap=ncmap,norm=norm1)

    axe.set_xlim([ilon[0,0],ilon[-1,-1]])
    axe.set_ylim([ilat[0,0],ilat[-1,-1]])
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

    axe.set_xlim([ilon[0,0],ilon[-1,-1]])
    axe.set_ylim([ilat[0,0],ilat[-1,-1]])
    axe.set_xticks(np.arange(np.ceil(ilon[0,0]),np.ceil(ilon[-1,-1]),0.5), crs=ccrs.PlateCarree())
    axe.set_yticks(np.arange(np.ceil(ilat[0,0]),np.ceil(ilat[-1,-1]),0.5), crs=ccrs.PlateCarree())
    axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
    axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))
    plt.savefig(figname,bbox_inches="tight",pad_inches=0.)
    del fig, axe

if __name__=='__main__':
    main_run()
