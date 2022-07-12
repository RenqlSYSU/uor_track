import xarray as xr
import numpy as np
import subprocess
import os 
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
from cartopy.io.shapereader import Reader
import cmaps

title_font=14
label_font=10
plt.rcParams["font.weight"] = "bold"
font = {'family': 'sans-serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black',
        }

def main_run():
    time = ['2021091500','2021102100']
    figdir = "/home/lzhenn/cooperate/fig/"
    '''
    for t in time:
        hycom = '/home/lzhenn/drv_field/hycom_subset/%s/hycom_glby_930_%s.nc'%(t,t)
        ds_hycom = xr.open_dataset(hycom)
        var = ds_hycom['surf_el'][0,:,:]
        test_map('%s/%s.jpg'%(figdir,t), '%s zeta'%t, [-1.6,1.6,0.2], 
            cmaps.BlueDarkRed18, var.data, var.lat, var.lon)
    '''
    roms_grid = '/home/lzhenn/Njord/Projects/Njord/roms_swan_grid/roms_d01_lp0d1.nc'
    ds_grid  = xr.open_dataset(roms_grid) # contain horizontal grid information
    test_map('%s/roms_h.jpg'%(figdir), 'bathymetry', [10,2530,20], 
        cmaps.MPL_rainbow, ds_grid['h'],
        ds_grid['lat_rho'],ds_grid['lon_rho'])

def test_map(figname, title, cnlev, fcolors, var,ilat,ilon):
    print(var.min())
    print(var.max())
    fig = plt.figure(figsize=(12,9),dpi=300)
    axe = plt.axes(projection=ccrs.PlateCarree())
    axe.set_title(title, fontsize=title_font, fontdict=font)
    #coast_shp2 = Reader(os.getenv("SHP_LIB")+
    #    "/china_coast/china_coastline.dbf").geometries()
    #coastline2 = cfeat.ShapelyFeature(coast_shp2, 
    #    ccrs.PlateCarree(), edgecolor="black",facecolor="none")
    #axe.add_feature(coastline2, linewidth=0.8,zorder=2)
    axe.add_feature(cfeat.NaturalEarthFeature('physical', 'land', '50m',
        edgecolor='k',facecolor=cfeat.COLORS['land']),linewidth=0.8)

    cnlevels = np.arange(cnlev[0], cnlev[1], cnlev[2])
    norm1 = colors.BoundaryNorm(boundaries=cnlevels, 
        ncolors=fcolors.N,extend='both')
    
    cont = axe.contourf(ilon, ilat, var, cnlevels,
         transform=ccrs.PlateCarree(),cmap=fcolors,extend="both",norm=norm1)

    axe.set_xlim([ilon.min(),ilon.max()])
    axe.set_ylim([ilat.min(),ilat.max()])
    axe.set_xticks(np.arange(np.ceil(ilon.min()),np.ceil(ilon.max()),5), crs=ccrs.PlateCarree())
    axe.set_yticks(np.arange(np.ceil(ilat.min()),np.ceil(ilat.max()),5), crs=ccrs.PlateCarree())
    axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
    axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))
    
    plt.colorbar(cont, ax=axe, shrink=.8)
    plt.savefig(figname,bbox_inches="tight",pad_inches=0.)
    del fig, axe

if __name__=='__main__':
    main_run()

