#!/usr/bin/env python
import cf
import cfplot as cfp
import xarray as xr
import numpy as np
import sys, os, subprocess, linecache, gc
import cartopy.crs as ccrs

def calc_month( frae, month, nday, filname):
    var  = np.empty( [len(month)], dtype=int ) # 4 or 12 
    datafile = "/home/users/qd201969/uor_track/fig/match_ens1_yes.dat"
    for nm in range(0,len(month),1):
        fras = frae+1
        frae = fras+24*nday[nm]-1
        print("month=%s, frame_s=%d, frame_e=%d" %(month[nm],fras,frae))
        rf = open("/home/users/qd201969/uor_track/fig/region.dat","w")
        rf.write("0 360\n-90 90\n"+str(fras)+" "+str(frae))
        rf.close()

        ret=subprocess.Popen("/home/users/qd201969/TRACK-1.5.2/utils/bin/censemble2 "+\
                filname+" "+filname+" 0 100 10 1 0 0 0 0 s 0 1 > rec",shell=True)
        ret.wait()

        a=linecache.getline(datafile, 4).strip().split(" ",1)
        term = a[1].strip().split(" ",1)
        var[nm] = int(term[0])
        linecache.clearcache()
    return var

def draw_month_traj(frae, month, nday, filname):
    var  = np.empty( [len(month)], dtype=int ) # 4 or 12
    fterm = filname.split("_")
    datafile = "/home/users/qd201969/uor_track/fig/match_ens1_yes.dat"
    figname = "traj-mon_"+fterm[1]+"_"+fterm[4]+"_"+fterm[-1]+".png"
    
    f0=cf.read("/home/users/qd201969/gtopo30_0.9x1.25.nc")
    phis=f0[2]
    phis=phis/9.8 # transfer from m2/s2 to m

    mlonl=0  #0  #
    mlonr=150#360#
    mlats=15 #0  #
    mlatn=70 #90 #
    ds = xr.open_dataset('/home/users/qd201969/data/ERA5_mon_u_1979-2020.nc')
    lat = ds.latitude
    lon = ds.longitude
    ilon = lon[(lon>=mlonl) & (lon<=mlonr)]
    ilat = lat[(lat>=mlats) & (lat<=mlatn)]
    da = ds['u'].sel(level=200,longitude=ilon,latitude=ilat).load()
    # increased performance by loading data into memory first, e.g., with load()
    uwnd = da.groupby(da.time.dt.month).mean('time')
    del ds, da
    gc.collect()
   
    cfp.setvars(file=figname)
    cfp.gopen(figsize=[20, 20],rows=4,columns=3,wspace=0.1,hspace=0.015,bottom=0.5)
    cfp.mapset(lonmin=mlonl, lonmax=mlonr, latmin=mlats, latmax=mlatn)
    for nm in range(0,len(month),1):
        fras = frae+1
        frae = fras+24*nday[nm]-1
        rf = open("/home/users/qd201969/uor_track/fig/region.dat","w")
        rf.write("0 360\n-90 90\n"+str(fras)+" "+str(frae))
        rf.close()
        
        ret=subprocess.Popen("/home/users/qd201969/TRACK-1.5.2/utils/bin/censemble2 "+\
                filname+" "+filname+" 0 100 10 1 0 0 0 0 s 0 1 > rec",shell=True)
        ret.wait()
        a=linecache.getline(datafile, 4).strip().split(" ",1)
        term = a[1].strip().split(" ",1)
        var[nm] = int(term[0])
        linecache.clearcache()
        print("month=%s, frame_s=%d, frame_e=%d, number= %d" %(month[nm],fras,frae,var[nm]))
        
        cfp.gpos(nm+1)
        cfp.levs(manual=[1500,3000,4500])
        cfp.con(phis,fill=False, lines=True, colors='k',linewidths=2,line_labels=False,
                title=fterm[-1]+" "+fterm[1]+"hPa "+fterm[4]+" "+month[nm]+" "+str(var[nm]))
        cfp.levs(manual=[30,32]) # jet stream
        cfp.con(uwnd[nm,:,:],x=ilon, y=ilat, ptype=1, fill=False,line_labels=False,
                lines=True, colors='indigo',linewidths=3.5)
        for nrr in range(0,len(lats),1):
            cfp.plotvars.mymap.plot([lonl[nrr],lonl[nrr],lonr[nrr],lonr[nrr],lonl[nrr]],
                    [latn[nrr],lats[nrr],lats[nrr],latn[nrr],latn[nrr]], 
                    linewidth=4, color='k', transform=ccrs.PlateCarree()) # filter box
        
        if var[nm] == 0:
            continue
        
        ret=subprocess.Popen("/home/users/qd201969/TRACK-1.5.2/utils/bin/tr2nc \
                "+datafile+" s /home/users/qd201969/TRACK-1.5.2/utils/TR2NC/tr2nc.meta",shell=True)
        ret.wait()
        f=cf.read(datafile+'.nc')
        g = f[2]
        g = g*1e5
        cfp.cscale('precip2_17lev')
        cfp.levs(min=0.0, max=8.0, step=0.5)
        cfp.traj(g, zorder=0, legend_lines=True, colorbar=False, linewidth=1.5)
        
    cfp.cscale('precip2_17lev')
    cfp.levs(min=0.0, max=8.0, step=0.5)
    cfp.cbar(position=[0.2, 0.48, 0.6, 0.01], title='Relative Vorticity (Hz)*1e5')
    cfp.gclose()
    subprocess.run("mogrify -bordercolor white -trim ./"+figname,shell=True) 
    return var

