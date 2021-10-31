'''
read  

20211026
'''

import xarray as xr
import numpy as np
import subprocess
import os 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import colors
from matplotlib.colors import ListedColormap, LightSource
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
from cartopy.io.shapereader import Reader
import cmaps

def sph2cart(r, theta, phi):
    '''spherical to cartesian transformation.'''
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z

def sphview(ax):
    '''returns the camera position for 3D axes in spherical coordinates'''
    r = np.square(np.max([ax.get_xlim(), ax.get_ylim()], 1)).sum()
    theta, phi = np.radians((90-ax.elev, ax.azim))
    return r, theta, phi

def ravzip(*itr):
    '''flatten and zip arrays'''
    return zip(*map(np.ravel, itr))

def draw_3d(title,lon,lat,var,norm,ncmap):
    fig = plt.figure(figsize=(12,9),dpi=300)
    axe = plt.axes(projection="3d")
    axe.set_title(title,fontsize=MIDFONT) 
    
    #ls = LightSource(270, 45)
    #rgb = ls.shade(var, cmap=ncmap, norm=norm, vert_exag=0.1, blend_mode='soft')
    #surf = axe.plot_surface(lon,lat,var,facecolors=rgb,rstride=1,cstride=1, alpha=1.,
    #                        linewidth=0, antialiased=False)
    #surf = axe.plot_surface(lon,lat,var,facecolor="lightgray",rstride=1,cstride=1, alpha=1.,
    #                        linewidth=0, antialiased=False)
    #surf = axe.plot_surface(lon,lat,var,cmap=ncmap,norm=norm,rstride=1,cstride=1, alpha=1.,
    #                        linewidth=0, antialiased=False)
    #plt.colorbar(surf, ax=axe, shrink=.8)
    
    #axe.set_facecolor('k')
    axe.set_zlim(-80,620)
    #axe.set_zscale('log')
    axe.view_init(elev=53, azim=-75)
    #axe.grid(False)
    #plt.axis('off')
    
    dx = np.diff(lon[1,1:3])
    dy = np.diff(lat[1:3,1])
    res = len(lon[0,:])
    xyz = np.array(sph2cart(*sphview(axe)), ndmin=3).T   #camera position in xyz
    zo = np.multiply([lon, lat, np.zeros_like(var)], xyz).sum(0)  #"distance" of bars from camera
    bars = np.empty(lon.shape, dtype=object)
    for i,(x,y,dz,o) in enumerate(ravzip(lon, lat, var, zo)):
        for nll in range(0,len(cnlevels),1):
            if nll==0 and dz < cnlevels[nll]:
                color0 = ncmap([nll])
                break
            elif nll<(len(cnlevels)-1) and dz >= cnlevels[nll] and dz < cnlevels[nll+1]:
                color0 = ncmap([nll+1])
                break
            else:
                color0 = ncmap([-1])

        j, k = divmod(i, res)
        bars[j, k] = pl = axe.bar3d(x, y, 0, dx, dy, dz, color0)
        pl._sort_zpos = o
   
    colourMap = plt.cm.ScalarMappable(cmap=ncmap,norm=norm)
    fig.colorbar(colourMap,extend='both',shrink=0.6)#,label='m'
    plt.savefig("%s%s_3d.jpg"%(figdir,title),bbox_inches="tight",pad_inches=0.1)
    del fig, axe

def draw_one(title,lon,lat,var,norm,ncmap):
    fig = plt.figure(figsize=(12,9),dpi=300)
    axe = plt.axes(projection=ccrs.PlateCarree())
    axe.set_title(title,fontsize=MIDFONT) 
    
    coast_shp = Reader(os.getenv("SHP_LIB")+"/china_coast/china_coastline.dbf").geometries()
    coastline = cfeat.ShapelyFeature(coast_shp, ccrs.PlateCarree(), edgecolor="black", facecolor="none")
    axe.add_feature(coastline, linewidth=0.8,zorder=1)

    #cont = axe.pcolormesh(lon, lat, var, 
    #        transform=ccrs.PlateCarree(),cmap=cmaps.precip2_17lev,norm=norm)
    cont = axe.pcolormesh(lon, lat, var, 
            transform=ccrs.PlateCarree(),cmap=ncmap,norm=norm)

    # Set the map bounds
    axe.set_xticks(np.arange(round(lon[0,0],1),round(lon[-1,-1],1),lon_sp), crs=ccrs.PlateCarree())
    axe.set_yticks(np.arange(round(lat[0,0],1),round(lat[-1,-1],1),lat_sp), crs=ccrs.PlateCarree())
    axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=""))
    axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=""))

    plt.colorbar(cont, ax=axe, shrink=.8)
    plt.savefig("%s%s.jpg"%(figdir,title),bbox_inches="tight",pad_inches=0.1)
    del fig, axe

BIGFONT=22
MIDFONT=18
SMFONT=14
lat_sp = 0.1 #1.0 #
lon_sp = 0.1 #1.0 #
figdir = "/home/lzhenn/cooperate/fig/"

filname=["/home/lzhenn/cooperate/data/Whole_HK_DTM_100m.nc",\
         "/home/lzhenn/njord_implement/domaindb/gba_norm/roms_d01.nc",\
         "/home/metctm1/array/data/Calypso/roms_d02.nc",\
         "/home/metctm1/array/data/Calypso/roms_d03.nc"]
fileout="/home/lzhenn/cooperate/data/dtm_bathymetry_d03.nc"
title1 = ["100m","2.2km","500m","100m"]

f  = xr.open_dataset(filname[3])
#var1 = f['h']
#var1 = f['mask_rho']
lats = f['lat_rho'].sel(eta_rho=0 ,xi_rho=1).data
latn = f['lat_rho'].sel(eta_rho=-1,xi_rho=1).data
lonl = f['lon_rho'].sel(eta_rho=1,xi_rho=0).data
lonr = f['lon_rho'].sel(eta_rho=1,xi_rho=-1).data

'''
print(var1.name,var1.sizes)
print(var0.name,var0.sizes)
var = var1
var[np.where(np.array([var1==1, var0==0]).any(axis=0))] = 0 # water
#var[np.where(var1==1 or var0==0] = 0 # water
'''

ncmap = ListedColormap(cmaps.GMT_globe(range(0,209,2))) #cmaps.precip2_17lev
#ncmap = ListedColormap(["red","blue","green"])
cnlevels = np.concatenate((np.arange(-63,0,1),np.arange(0,615,15)))
norm  = colors.BoundaryNorm(boundaries=cnlevels, ncolors=ncmap.N,extend="both")
#norm = colors.TwoSlopeNorm(vmin=-60, vcenter=0, vmax=450)
#draw_one("100m_terrain-bathy",lon,lat,var,norm,ncmap)
#draw_3d("100m_terrain-bathy",lon,lat,var,norm,ncmap)

for nf in range(1,4,1):#len(filname),1):
    f  = xr.open_dataset(filname[nf])
    lat = f['lat_rho'].sel(xi_rho=1).data
    lon = f['lon_rho'].sel(eta_rho=1).data
    ixc = f.xi_rho[(lon>=lonl) & (lon<=lonr)]
    iyc = f.eta_rho[(lat>=lats) & (lat<=latn)]
    del lat, lon
    var1 = f['h'].sel(eta_rho=iyc, xi_rho=ixc).data
    ilat = f['lat_rho'].sel(eta_rho=iyc, xi_rho=ixc).data
    ilon = f['lon_rho'].sel(eta_rho=iyc, xi_rho=ixc).data

    f  = xr.open_dataset(filname[0])
    var0 = f['dtm'].sel(lat=ilat[:,1],lon=ilon[1,:],method="nearest")
    var = np.subtract(var0.data,var1.data)

    print("%s resolution: %f, %f"%(filname[nf],(ilat[0,0]-ilat[1,1]),(ilon[0,0]-ilon[1,1])))
    #term = filname[nf].split("/")
    #title=term[-1].split(".")
    
    #draw_one(title[0],lon,lat,var,norm,ncmap)
    #draw_one("%s_terrain-bathy"%title1[nf],ilon,ilat,var,norm,ncmap)
    draw_3d("%s_terrain-bathy"%title1[nf],ilon,ilat,var,norm,ncmap)

