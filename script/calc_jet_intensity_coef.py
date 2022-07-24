import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

title_font=18
label_font=14
plt.rcParams["font.weight"] = "bold"
font = {'family': 'sans-serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black', 
        }

path = '/WORK/sysu_hjkx_ys/renql/JET_NUDG/input2nudg'
case = ['CTRL','NUDG']
ftime  = pd.date_range(start='2019-01-01 00',end='2019-12-31 23',freq='1D')
print(ftime)

def main_run():
    #for nm in [7,8]:
    #    draw_map()
    #    draw_vertical(nm)
    #draw_vertical(6)
    #draw_lat_time(200)
    draw_map(6,200)
    
def calc_jet_inte():
    lonl = 105
    lonr = 130
    lats = 20
    latn = 50

    ds1 = xr.open_dataset('%s/F2000_CAM5.cam.h1.ESM.clim.U.nc'%path)
    lat = ds1['lat'].data
    lon = ds1['lon'].data
    ilon = lon[(lon>=lonl) & (lon<=lonr)]
    ilat = lat[(lat>=lats) & (lat<=latn)]
    da1 = ds1['U'].sel(lon=ilon,lat=ilat).mean('lon').max(['lev','lat'])
    print(da1)
    print("")

    ds2 = xr.open_dataset('%s/F2000_CAM5.cam.h1.ESM.clim.month_U.nc'%path)
    da2 = ds2['U'].sel(lon=ilon,lat=ilat).mean('lon').max()
    print(da2)
    print("")

    coe = da1/da2
    print(coe)
    print("")

    draw_ts(coe)

    var = ds1['U']
    var.data = (coe*ds2['U']).data
    print(var)
    print("")
    ds3 = var.to_dataset(name="U")
    ds3.attrs["description"] = 'fixed jet with annual variation intensity'
    ds3.to_netcdf('%s/F2000_CAM5.cam.h1.ESM.clim.daily_U.nc'%path,"w")

def draw_lat_time(lev):
    lonl = 105
    lonr = 130
    lats = 20
    latn = 50
    ds1 = xr.open_dataset('%s/F2000_CAM5.cam.h1.ESM.clim.U.nc'%path)
    lat = ds1['lat'].data
    lon = ds1['lon'].data
    ilon = lon[(lon>=lonl) & (lon<=lonr)]
    ilat = lat[(lat>=lats) & (lat<=latn)]
    da1 = ds1['U'].sel(lev=lev,lon=ilon,lat=ilat,method='nearest').mean('lon')
    print(da1)
    print("")

    ds2 = xr.open_dataset('%s/F2000_CAM5.cam.h1.ESM.clim.daily_U.nc'%path)
    da2 = ds2['U'].sel(lev=lev,lon=ilon,lat=ilat,method='nearest').mean('lon')
    print(da2)
    var = [da1.data,da2.data]
    del ds1, ds2

    cnlevels = np.arange(-3,54,3)
    fcolors = plt.cm.get_cmap('rainbow',18) 
    norm = colors.BoundaryNorm(boundaries=cnlevels, 
        ncolors=fcolors.N)

    fig = plt.figure(figsize=(12,12),dpi=300)
    ax = fig.subplots(2,1)
    for nr in range(0,2,1):
        axe = ax[nr]
        axe.set_title('%s %d %d-%dE'%(case[nr],lev,lonl,lonr),fontsize=title_font, fontdict=font)
        cont = axe.contourf(ftime, da1.lat, var[nr].transpose(), 
            cnlevels, cmap=fcolors, norm=norm)
        axe.set_ylabel("Lat",fontsize=title_font, fontdict=font)  # Add a y-label to the axes.
        axe.set_xlabel("time",fontsize=title_font, fontdict=font)  # Add a y-label to the axes.
        axe.tick_params(axis='both', which='major', labelsize=title_font)
        axe.xaxis.set_major_formatter(mdates.DateFormatter("%m"))
        plt.colorbar(cont, ax=axe, orientation='horizontal')
    plt.savefig('%s/jet_lat_time_%d.jpg'%(path,lev), bbox_inches='tight',pad_inches=0.01)

def draw_map(nm,lev):
    lonl = 70
    lonr = 160
    lats = 10
    latn = 60
    lat_sp = 10
    lon_sp = 20 #60 #
    ds1 = xr.open_dataset('%s/F2000_CAM5.cam.h1.ESM.clim.U.nc'%path)
    lat = ds1['lat'].data
    lon = ds1['lon'].data
    ilon = lon[(lon>=lonl) & (lon<=lonr)]
    ilat = lat[(lat>=lats) & (lat<=latn)]
    da1 = ds1['U'].sel(time=ds1.time.dt.month.isin(nm),lev=lev,
        lon=ilon,lat=ilat).mean(['time'])
    print(da1)
    print("")

    ds2 = xr.open_dataset('%s/F2000_CAM5.cam.h1.ESM.clim.daily_U.nc'%path)
    da2 = ds2['U'].sel(time=ds2.time.dt.month.isin(nm),lev=lev,
        lon=ilon,lat=ilat).mean(['time'])
    #ds2 = xr.open_dataset('%s/F2000_CAM5.cam.h1.ESM.clim.month_U.nc'%path)
    #da2 = ds2['U'].sel(lon=ilon,lat=ilat).mean('lon')
    print(da2)
    var = [da1.data,da2.data]
    del ds1, ds2

    cnlevels = np.arange(-3,54,3)
    fcolors = plt.cm.get_cmap('rainbow',18) 
    norm = colors.BoundaryNorm(boundaries=cnlevels, 
        ncolors=fcolors.N)
    
    topofile = '/WORK/sysu_hjkx_ys/cesm/input/atm/cam/topo/USGS-gtopo30_0.9x1.25_remap_c051027.nc'
    ds = xr.open_dataset(topofile)
    phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").load()
    phis = phis/9.8 # transfer from m2/s2 to m
    del ds
    
    fig = plt.figure(figsize=(12,12),dpi=300)
    ax = fig.subplots(1,2)
    for nr in range(0,2,1):
        axe = ax[nr]
        axe.set_title('%s %d %d'%(case[nr],nm,lev),
            fontsize=title_font, fontdict=font)
        axe.add_feature(cfeat.GSHHSFeature(levels=[1,2],edgecolor='k')
                , linewidth=0.8, zorder=1)
        cont = axe.contourf(ilon, ilat, var[nr], cnlevels, 
             transform=ccrs.PlateCarree(), cmap=fcolors, norm=norm)
        topo = axe.contour(ilon, ilat, phis, [1500,3000], 
             transform=ccrs.PlateCarree(),colors='black',linewidths=1.5)
        axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
        axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
        axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
        axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))
        plt.colorbar(cont, ax=axe, orientation='horizontal')
    plt.savefig('%s/jet_%d.jpg'%(path,nm), bbox_inches='tight',pad_inches=0.01)
    
def draw_vertical(nm):
    lonl = 105
    lonr = 130
    lats = 20
    latn = 50
    ds1 = xr.open_dataset('%s/F2000_CAM5.cam.h1.ESM.clim.U.nc'%path)
    lat = ds1['lat'].data
    lon = ds1['lon'].data
    ilon = lon[(lon>=lonl) & (lon<=lonr)]
    ilat = lat[(lat>=lats) & (lat<=latn)]
    da1 = ds1['U'].sel(time=ds1.time.dt.month.isin(nm),
        lon=ilon,lat=ilat).mean(['lon','time'])
    print(da1)
    print("")

    #ds2 = xr.open_dataset('%s/F2000_CAM5.cam.h1.ESM.clim.daily_U.nc'%path)
    #da2 = ds2['U'].sel(time=ds2.time.dt.month.isin(nm),
    #    lon=ilon,lat=ilat).mean(['lon','time'])
    ds2 = xr.open_dataset('%s/F2000_CAM5.cam.h1.ESM.clim.month_U.nc'%path)
    da2 = ds2['U'].sel(lon=ilon,lat=ilat).mean('lon')
    print(da2)
    var = [da1.data,da2.data]
    del ds1, ds2

    cnlevels = np.arange(-3,54,3)
    fcolors = plt.cm.get_cmap('rainbow',18) 
    norm = colors.BoundaryNorm(boundaries=cnlevels, 
        ncolors=fcolors.N)

    fig = plt.figure(figsize=(12,12),dpi=300)
    ax = fig.subplots(1,2)
    for nr in range(0,2,1):
        axe = ax[nr]
        axe.set_title('%s %d %d-%dE'%(case[nr],nm,lonl,lonr),
            fontsize=title_font, fontdict=font)
        cont = axe.contourf(da1.lat, da1.lev, var[nr], cnlevels, 
             cmap=fcolors, norm=norm)
        axe.set_xlabel('Lat',fontsize=label_font,fontdict=font)
        axe.set_ylabel('sigma lev',fontsize=label_font,fontdict=font)
        axe.set_yscale('symlog')
        axe.set_yticklabels(np.arange(1000,5,-100))
        axe.set_yticks(np.arange(1000,5,-100))
        axe.set_ylim(da1.lev.max(),100)
        plt.colorbar(cont, ax=axe, orientation='horizontal')
    
    plt.savefig('%s/jet_%d.jpg'%(path,nm), bbox_inches='tight',pad_inches=0.01)

def draw_ts(ts):
    fig = plt.figure(figsize=(9,9),dpi=200)
    ax  = fig.subplots(1, 1) #sharex=True, sharey=True
    ax.set_title("coe" ,fontsize=title_font,fontdict=font)

    ax.plot(ftime,ts,linewidth=2)
    ax.set_xlabel("time",fontsize=title_font)  # Add a y-label to the axes.
    ax.tick_params(axis='both', which='major', labelsize=title_font)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))

    fig.tight_layout(w_pad=0.5,h_pad=1) #,rect=(0,bmlo,1,1)
    fig.savefig("%s/coe_annual_ts.png"%(path))

if __name__=='__main__':
    main_run()
