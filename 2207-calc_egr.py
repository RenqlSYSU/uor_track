#!/usr/bin/env python
import sys
import subprocess
import xarray as xr
import numpy as np
import gc #garbage collector
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.lines import Line2D
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmaps
#matplotlib.use('Agg')

title_font=14
label_font=10
plt.rcParams["font.weight"] = "bold"
font = {'family': 'sans-serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black', 
        }

level = [225, 250, 300, 450, 500, 550, 825, 850, 875]
path = '/gws/nopw/j04/ncas_generic/users/renql/ERA5_subdaily'
outdir = '/home/users/qd201969/uor_track/mdata'
figdir = '/home/users/qd201969/uor_track/fig'
lonl=15 
lonr=145
lats=15
latn=70
lat_sp = 20
lon_sp = 30 #60 #
ga = 9.80665 # Gravitational acceleration

def main_run():
    #calc_monthly_egr()
    draw_season_4x3()

def calc_monthly_egr():
    ds = xr.open_dataset("%s/u/ERA5_u_1980.nc"%path)
    var = ds['u'].sel(level=level).load()
    var = var.groupby(var.time.dt.month).mean('time')
    var.data = np.zeros(var.data.shape)

    lat = ds['latitude'].data
    del ds
    f0  = 2*(2*np.pi/24.0/3600.0)*np.sin(lat*np.pi/180.0)
    nyear = 0
    for year in range(1980,2000,1):
        print('calc for %d'%year)
 #       ds = xr.open_dataset("%s/t/ERA5_t_%d.nc"%(path,year))
 #       t = ds['t'].sel(level=level).groupby(var.time.dt.month).mean('time').data
 #       ds = xr.open_dataset("%s/z/ERA5_z_%d.nc"%(path,year))
 #       z = ds['z'].sel(level=level).groupby(var.time.dt.month).mean('time').data/ga
 #       ds = xr.open_dataset("%s/u/ERA5_u_%d.nc"%(path,year))
 #       u = ds['u'].sel(level=level).groupby(var.time.dt.month).mean('time')
        ds = xr.open_dataset("%s/t/ERA5_t_%d.nc"%(path,year))
        t = ds['t'].sel(level=level).data
        ds = xr.open_dataset("%s/z/ERA5_z_%d.nc"%(path,year))
        z = ds['z'].sel(level=level).data/ga
        ds = xr.open_dataset("%s/u/ERA5_u_%d.nc"%(path,year))
        u = ds['u'].sel(level=level)
        del ds
        u.data = calc_egr(t, z, u.data, f0)
        u = u.groupby(u.time.dt.month).mean('time')
        var.data = var.data + u.data 
        nyear = nyear + 1

    print('total year %d'%nyear)
    var.data = var.data*86400/nyear
    var.attrs['long_name'] = 'Eady growth rate'
    var.attrs['units'] = '1/day'
    print(var)
    ds1 = var.to_dataset(name="egr")
    ds1.to_netcdf('%s/month_egr_6h_0.25.nc'%outdir,'w')

def calc_egr(t,z,u,f0):
    t = calc_pot_temp(t, np.array(level), 1)
    brunt = np.sqrt(center_diff(t,z)*ga/t) 
    # negative is nan

    brunt = np.ma.array(brunt,mask=(brunt<1e-10))
    f00 = np.broadcast_to(f0.reshape(len(f0), 1), u.shape)
    term = 0.3098*np.abs(center_diff(u,z))*f00/brunt
    return term

def calc_pot_temp(t,p,dim):
    # dim indicating which dimension of t 
    # is similar to p
    term = np.moveaxis(t,dim,-1)*np.power(1000/p,0.286)
    return np.moveaxis(term,-1,dim)

def center_diff(var,z):
    term = var
    term[:,1:-1,:,:] = np.divide(var[:,0:-2,:,:]-var[:,2:,:,:],
        z[:,0:-2,:,:]-z[:,2:,:,:])
    term[:,0,:,:] = np.divide(var[:,0,:,:]-var[:,1,:,:],
        z[:,0,:,:]-z[:,1,:,:])
    term[:,-1,:,:] = np.divide(var[:,-2,:,:]-var[:,-1,:,:],
        z[:,-2,:,:]-z[:,-1,:,:])
    return term 

def draw_season_4x3():
    lev = [850,500,250]
    titls= ['DJF','MAM','JJA','SON']
    cnlev = np.arange(0,0.7,0.05)
    
    ds = xr.open_dataset('%s/month_egr_6h_0.25.nc'%outdir)
    lat = ds.latitude
    lon = ds.longitude
    ilon = lon[(lon>=lonl) & (lon<=lonr)]
    ilat = lat[(lat>=lats) & (lat<=latn)]
    var = 100*ds['egr'].sel(level=lev,longitude=ilon,latitude=ilat).data
    
    nrow = 4 #6 #
    ncol = 3 #2 #
    bmlo = 0.4 #0.25 #
    
    fig = plt.figure(figsize=(12,12),dpi=300)
    ax = fig.subplots(nrow, ncol, subplot_kw=dict(
        projection=ccrs.PlateCarree(central_longitude=180.0))) #sharex=True, sharey=True

    ncmap = colors.ListedColormap(cmaps.topo_15lev(range(0,16,1))[::-1])
    norm = colors.BoundaryNorm(boundaries=cnlev,
        ncolors=ncmap.N,extend='both')
    
    uwndpath = '/gws/nopw/j04/ncas_generic/users/renql/ERA5_mon/ERA5_mon_u_1979-2020.nc'
    ds = xr.open_dataset(uwndpath)
    da = ds['u'].sel(level=200,longitude=ilon,
        latitude=ilat,method="nearest").load()
    uwnd = da.groupby(da.time.dt.month).mean('time').data
    del ds, da

    ds = xr.open_dataset("/home/users/qd201969/gtopo30_0.9x1.25.nc")
    phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").data
    phis = phis/9.8 # transfer from m2/s2 to m
    del ds
    gc.collect()

    for nl in range(0,3,1):
        for nm in range(0,nrow,1):
            if nm == 0:
                shad = (var[0,nl,:,:]+var[1,nl,:,:]+var[11,nl,:,:])/3.0
                uwnd1 = (uwnd[0,:,:]+uwnd[1,:,:]+uwnd[11,:,:])/3.0
            else:
                shad = np.mean(var[(3*nm-1):(3*nm+2),nl,:,:],axis=0)
                uwnd1 = np.mean(uwnd[(3*nm-1):(3*nm+2),:,:],axis=0)
            print('%d EGR %s : min = %f ; max = %f'%(lev[nl],titls[nm],
                np.nanmin(shad),np.nanmax(shad)))
            axe = ax[nm][nl]
            axe.add_feature(cfeat.GSHHSFeature(levels=[1,2],
                edgecolor='k'), linewidth=0.8, zorder=1)
            axe.set_title("%dhPa %s EGR"%(lev[nl],titls[nm]
                ),fontsize=title_font,fontdict=font)

            cont = axe.contourf(ilon, ilat, shad, cnlev, 
                 transform=ccrs.PlateCarree(),cmap=ncmap,extend='both',norm=norm)
            topo = axe.contour(ilon, ilat, phis, [1500,3000], 
                 transform=ccrs.PlateCarree(),colors='black',linewidths=1.5)
            jets = axe.contour(ilon, ilat, uwnd1, [30,40,50], 
                 transform=ccrs.PlateCarree(),colors='r',linewidths=2.2)
            
            if nl == 0:
                axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
                axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
            if nm == (nrow-1):
                axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
                axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))

    position = fig.add_axes([0.45, bmlo+0.005, 0.45, 0.01]) #left, bottom, width, height
    cb = plt.colorbar(cont, cax=position ,orientation='horizontal')#, shrink=.9)
    plt.figtext(0.35,bmlo-0.005, 'EGR (1/day)',fontsize=title_font,
        horizontalalignment='left',verticalalignment='bottom')
    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig('%s/seasonal_egr.png'%(figdir), bbox_inches='tight',pad_inches=0.01)

if __name__=='__main__':
    main_run()

