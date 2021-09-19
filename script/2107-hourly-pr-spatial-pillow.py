import xarray as xr
import numpy as np
import subprocess
import os 
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
from cartopy.io.shapereader import Reader
import cmaps
from PIL import Image, ImageDraw, ImageSequence
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                         cartopy_ylim, latlon_coords)
from copy import copy
import shapely.geometry as sgeom

def find_side(ls, side):
    """
    Given a shapely LineString which is assumed to be rectangular, return the
    line corresponding to a given side of the rectangle.
    """
    minx, miny, maxx, maxy = ls.bounds
    points = {'left': [(minx, miny), (minx, maxy)],
            'right': [(maxx, miny), (maxx, maxy)],
            'bottom': [(minx, miny), (maxx, miny)],
            'top': [(minx, maxy), (maxx, maxy)],}
    return sgeom.LineString(points[side])

def lambert_xticks(ax, ticks):
    """Draw ticks on the bottom x-axis of a Lambert Conformal projection."""
    te = lambda xy: xy[0]
    lc = lambda t, n, b: np.vstack((np.zeros(n) + t, np.linspace(b[2], b[3], n))).T
    xticks, xticklabels = _lambert_ticks(ax, ticks, 'bottom', lc, te)
    ax.xaxis.tick_bottom()
    ax.set_xticks(xticks)
    ax.set_xticklabels([ax.xaxis.get_major_formatter()(xtick) for xtick in xticklabels])

def lambert_yticks(ax, ticks):
    """Draw ticks on the left y-axis of a Lamber Conformal projection."""
    te = lambda xy: xy[1]
    lc = lambda t, n, b: np.vstack((np.linspace(b[0], b[1], n), np.zeros(n) + t)).T
    yticks, yticklabels = _lambert_ticks(ax, ticks, 'left', lc, te)
    ax.yaxis.tick_left()
    ax.set_yticks(yticks)
    ax.set_yticklabels([ax.yaxis.get_major_formatter()(ytick) for ytick in yticklabels])

def _lambert_ticks(ax, ticks, tick_location, line_constructor, tick_extractor):
    """Get the tick locations and labels for an axis of a Lambert Conformal projection."""
    outline_patch = sgeom.LineString(ax.outline_patch.get_path().vertices.tolist())
    axis = find_side(outline_patch, tick_location)
    n_steps = 30
    extent = ax.get_extent(ccrs.PlateCarree())
    _ticks = []
    for t in ticks:
        xy = line_constructor(t, n_steps, extent)
        proj_xyz = ax.projection.transform_points(ccrs.Geodetic(), xy[:, 0], xy[:, 1])
        xyt = proj_xyz[..., :2]
        ls = sgeom.LineString(xyt.tolist())
        locs = axis.intersection(ls)
        if not locs:
            tick = [None]
        else:
            tick = tick_extractor(locs.xy)
        _ticks.append(tick[0])
    # Remove ticks that aren't visible:    
    ticklabels = copy(ticks)
    while True:
        try:
            index = _ticks.index(None)
        except ValueError:
            break
        _ticks.pop(index)
        ticklabels.pop(index)
    return _ticks, ticklabels

# constants
BIGFONT=22
MIDFONT=18
SMFONT=14

# define date of start time and forcast time
stime  = ['2021','06','22'] # year, month, date, hour
ftime1 = ['2021','06','23','08']
ftime2 = ['2021','06','24','08']

figtle = "Precipitation (mm/hour)"
case = ['Njord','PATH']
domin = ['03','04']
nc = 0

path = '/home/lzhenn/cooperate/data'
filedir= [path+'/Njord/'+''.join(stime)+'/',\
          path+'/PATH/'+''.join(stime[0:2])+'/'+''.join(stime)+'12/']
coast_shp = Reader(os.getenv('SHP_LIB')+'/china_coast/china_coastline.dbf').geometries()

# contour level for precip
cnlevels = [0.1, 0.5, 1, 2, 3, 5, 8, 12, 16, 20, 25, 30, 40, 50, 70, 100, 150]
norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=cmaps.precip2_17lev.N,extend='both')

fn_stream = subprocess.check_output('ls '+filedir[nc]+'wrfout_d'+domin[nc]+'_*', shell=True).decode('utf-8')
fn_list   = fn_stream.split()
print(fn_list[0])
print('filenumber : '+str(len(fn_list)))

ncfile = Dataset(fn_list[0])
var    = getvar(ncfile, "RAINC")
varnp1 = to_np(getvar(ncfile,"RAINC")) + to_np(getvar(ncfile,"RAINNC"))
nproj  = str(var.attrs['projection']).split('(')
print(nproj)

# Get the latitude and longitude points
lats, lons = latlon_coords(var)

# Get the cartopy mapping object
cart_proj = get_cartopy(var)

for itm in fn_list[1:]:
    fn_splt = itm.split('_')
    print(fn_splt) #['wrfout', 'd03', '2021-06-22', '11:00:00']
    figtle = case[nc]+' '+fn_splt[2]+' '+fn_splt[3]+' precip(mm/h)'
    figdir = '/home/lzhenn/cooperate/fig/'+case[nc]+'.S'+''.join(stime)+ \
             '.F'+''.join(fn_splt[2].split('-'))+'-'+fn_splt[3]+'.pr.png'

    # Open the NetCDF file
    ncfile = Dataset(itm)
    varnp  = to_np(getvar(ncfile,"RAINC")) + to_np(getvar(ncfile,"RAINNC"))-varnp1
    varnp1 = to_np(getvar(ncfile,"RAINC")) + to_np(getvar(ncfile,"RAINNC"))

    fig = plt.figure(figsize=(9,9))
        
    # Set the GeoAxes to the projection used by WRF
    axe = fig.add_axes([0.08, 0.05, 0.8, 0.94], projection=cart_proj)
    #axe = plt.axes(projection=cart_proj)
    #axe.add_geometries(coast_shp, ccrs.PlateCarree(),facecolor='none', edgecolor='grey',linewidth=.5, zorder = 1)
    #axe.add_feature(cfeat.COASTLINE.with_scale('10m'), linewidth=1,color='k')
    coast_shp = Reader(os.getenv('SHP_LIB')+'/china_coast/china_coastline.dbf').geometries()
    coastline = cfeat.ShapelyFeature(coast_shp, ccrs.PlateCarree(), edgecolor='grey', facecolor='none')
    axe.add_feature(coastline, linewidth=0.5)

    plt.contourf(to_np(lons), to_np(lats), varnp, cnlevels, 
                 transform=ccrs.PlateCarree(),cmap=cmaps.precip2_17lev,norm=norm,extend='both')

    plt.colorbar(ax=axe, shrink=.58)
    
    # Set the map bounds
    if nproj[0] == 'LambertConformal':
        fig.canvas.draw()
        # Define gridline locations and draw the lines using cartopy's built-in gridliner:
        xticks = list(np.arange(np.ceil(lons[0,0]+4.0),np.floor(lons[0,-1]+4.0),0.5))
        yticks = list(np.arange(np.ceil(lats[0,0]+4.0),np.floor(lats[-1,0]+4.0),0.5))
        axe.gridlines(xlocs=xticks, ylocs=yticks)
        # Label the end-points of the gridlines using the custom tick makers:
        axe.xaxis.set_major_formatter(LONGITUDE_FORMATTER) 
        axe.yaxis.set_major_formatter(LATITUDE_FORMATTER)
        # xticks and yticks only is list
        lambert_xticks(axe, xticks)  
        lambert_yticks(axe, yticks)
    else:
        axe.set_xlim(cartopy_xlim(var))
        axe.set_ylim(cartopy_ylim(var))
        axe.set_xticks(np.arange(np.ceil(lons[0,0]),np.floor(lons[0,-1]),1.0), crs=ccrs.PlateCarree())
        axe.set_yticks(np.arange(np.ceil(lats[0,0]),np.floor(lats[-1,0]),1.0), crs=ccrs.PlateCarree())
        axe.xaxis.set_major_formatter(LONGITUDE_FORMATTER) 
        axe.yaxis.set_major_formatter(LATITUDE_FORMATTER)
        axe.xaxis.set_major_formatter(LongitudeFormatter())
        axe.yaxis.set_major_formatter(LatitudeFormatter())
    
    plt.title(figtle,fontsize=SMFONT)
    plt.savefig(figdir)
    plt.show()
    del fig, axe

fn_stream = subprocess.check_output('ls /home/lzhenn/cooperate/fig/'+case[nc]+'.S'+''.join(stime)+'*.pr.png', shell=True).decode('utf-8')
fn_list   = fn_stream.split()
print(fn_list[0])
print('filenumber : '+str(len(fn_list)))
gif_name = './2.gif'

frames = []
for itm in fn_list:
    frame = Image.open(itm)
    frames.append(frame)

frames[0].save(gif_name, save_all=True, append_images=frames[1:],\
            duration = 1000, loop=0, disposal=1)
subprocess.run('rm -f /home/lzhenn/cooperate/fig/'+case[nc]+'.S'+''.join(stime)+'*.pr.png', shell=True)

