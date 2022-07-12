import xarray as xr
import pandas as pd
import numpy as np
import subprocess, os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import colors
from matplotlib import cm
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.io.shapereader import Reader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmaps, wrf
from netCDF4 import Dataset
from datetime import datetime

title_font=12
label_font=10
plt.rcParams["font.weight"] = "bold"
font = {'family': 'sans-serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black',
        }

figdir = "/home/lzhenn/cooperate/fig"
wrfpath = ['/home/lzhenn/cooperate/data/case_study/coupled',
           '/home/lzhenn/cooperate/data/case_study/wrfonly']
case = ['Coupled','Uncoupled']
slat = 21.56239
slon = 112.1174
elat = 22.62169
elon = 115.3909
startp = wrf.CoordPair(lat=slat,lon=slon)
endp = wrf.CoordPair(lat=elat,lon=elon)
angle = np.arctan((elat-slat)/(elon-slon)/np.cos(np.pi*22.0/180.0))
vlev = np.arange(100,1600,100) 

def main_run():
    fn_stream = subprocess.check_output('ls %s'%wrfpath[0],
        shell=True).decode('utf-8')
    stime0 = fn_stream.split()
    
    fn_stream = subprocess.check_output('ls %s'%wrfpath[1],
        shell=True).decode('utf-8')
    stime = fn_stream.split()
    #draw_ts(stime[0])
    draw_height_time(stime[0])
    '''
    for st in stime[1:]:
        if st in stime0:
            draw_ts(st)
            draw_height_time(st)
    '''
def draw_height_time(stime):
    title = 'onshore wind (m/s) %s'%stime
    figname = '%s/wind_%s.png'%(figdir,stime)
    cnlevels = np.arange(-6,6.1,0.5) #500Z, gpm
    #ncmap = colors.ListedColormap(cm.bwr(np.linspace(0, 1, 26)))#cmaps.BlueDarkRed18
    ncmap = colors.ListedColormap(np.vstack((cm.winter(np.linspace(0,1,13)),
        cm.autumn(np.linspace(0,1,13)[::-1]))))
    norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=ncmap.N, extend='both')
    
    fig = plt.figure(figsize=(12,9),dpi=300)
    #fig.suptitle(title,fontsize=title_font, fontdict=font)
    ax  = fig.subplots(2,1)
    for nc in range(len(case)):
        fn_list = subprocess.check_output('ls %s/%s/wrfout_d03_*'%(
            wrfpath[nc],stime), shell=True).decode('utf-8').split()
        print('%s filenumber=%d'%(fn_list[0],len(fn_list)))
        var = np.empty((len(fn_list),len(vlev)), dtype = float)
        ftime = []
        for nt in range(len(fn_list)):
            fn_splt = fn_list[nt].split('_')
            ftime.append(datetime.strptime('%s %s'%(fn_splt[-2],fn_splt[-1]),
                '%Y-%m-%d %H:00:00'))
            ncfile = Dataset(fn_list[nt])
            height = wrf.getvar(ncfile,'height')
            term = wrf.getvar(ncfile,'ua')
            u = wrf.vertcross(term, height, vlev, wrfin=ncfile, stagger='m', 
                start_point=startp, end_point=endp, latlon=True)
            term = wrf.getvar(ncfile,'va')
            v = wrf.vertcross(term, height, vlev, wrfin=ncfile, stagger='m', 
                start_point=startp, end_point=endp, latlon=True)
            var[nt,:] = np.mean(-u.data*np.sin(angle)+v.data*np.cos(angle),
                axis=1)
        print('min %f ; max %f'%(var.min(),var.max()))

        axe = ax[nc]
        axe.set_title("%s %s"%(case[nc],title), fontsize=title_font, fontdict=font)
        cont = axe.contourf(ftime, vlev, var.transpose(), cnlevels, 
             cmap=ncmap,extend='both',norm=norm)
        axe.set_ylabel("height (m)",fontsize=label_font, fontdict=font)  # Add an x-label to the axes.
        axe.tick_params(axis='both', which='major', labelsize=label_font)
        if nc==0:
            axe.set_xticks([])
        else:
            axe.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d_%H:%M"))
            plt.setp(axe.get_xticklabels(), rotation=15, ha="right")
    
    bmlo = 0.5
    position = fig.add_axes([0.2, bmlo+0.001, 0.7, 0.01]) #left, bottom, width, height
    cb = plt.colorbar(cont, cax=position ,orientation='horizontal')#, shrink=.9)
    #cb.set_label(label=var.long_name, size=title_font) #, weight='bold'
    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig(figname,bbox_inches="tight",pad_inches=0.)

def draw_ts(stime):
    title = '10m wind (m/s) %s'%stime
    figname = '%s/10mwind_%s.png'%(figdir,stime)
    
    fig = plt.figure(figsize=(12,9),dpi=300)
    axe = fig.add_axes([0.05, 0.05, 0.85, 0.45])
    axe.set_title(title, fontsize=title_font, fontdict=font)
    
    for nc in range(len(case)):
        fn_list = subprocess.check_output('ls %s/%s/wrfout_d03_*'%(
            wrfpath[nc],stime), shell=True).decode('utf-8').split()
        print('%s filenumber=%d'%(fn_list[0],len(fn_list)))
        var = np.empty((len(fn_list)), dtype = float)
        ftime = []
        for nt in range(len(fn_list)):
            fn_splt = fn_list[nt].split('_')
            ftime.append(datetime.strptime('%s %s'%(fn_splt[-2],fn_splt[-1]),
                '%Y-%m-%d %H:00:00'))

            ncfile = Dataset(fn_list[nt])
            term = wrf.getvar(ncfile,'U10')
            u = wrf.interpline(term, wrfin=ncfile, stagger='m', 
                start_point=startp, end_point=endp, latlon=True)
            term = wrf.getvar(ncfile,'V10')
            v = wrf.interpline(term, wrfin=ncfile, stagger='m', 
                start_point=startp, end_point=endp, latlon=True)
            var[nt] = np.mean(-u.data*np.sin(angle)+v.data*np.cos(angle))

        axe.plot(ftime, var, label=case[nc], linewidth=2)
    axe.tick_params(axis='both', which='major',labelsize=label_font) 
    axe.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d_%H:00"))
    plt.setp(axe.get_xticklabels(), rotation=30, ha="right")
    axe.grid(True, which="both", color='grey', linestyle='--', linewidth=1)
    axe.legend(fontsize=label_font)  # Add a legend.
    plt.savefig(figname,bbox_inches="tight",pad_inches=0.)

if __name__=='__main__':
    main_run()

