import xarray as xr
import numpy as np
import subprocess
import gc #garbage collector
import cf
import cfplot as cfp

lonl=0  #0  #
lonr=150#360#
lats=15 #0  #
latn=70 #90 #

path = '/home/users/qd201969/data/'
var_name = ['u','sp']
nv = 0

# -------------- read data and calc climatology ---------------
ds = xr.open_dataset(path+'ERA5_mon_'+var_name[nv]+'_1979-2020.nc')
lat = ds.latitude
lon = ds.longitude
ilon = lon[(lon>=lonl) & (lon<=lonr)]
ilat = lat[(lat>=lats) & (lat<=latn)]
if nv == 1:
    da = ds[var_name[nv]].sel(longitude=ilon,latitude=ilat).load()
    da.values /= 100.0 # convert Pa to hPa
    lev  =[400,1200,50]
else:
    da = ds[var_name[nv]].sel(level=200,expver=1,longitude=ilon,latitude=ilat).load()
    lev  =[20,52,2]
# increased performance by loading data into memory first, e.g., with load()

var = da.groupby(da.time.dt.month).mean('time')
print(var)
del ds
gc.collect()

f0=cf.read("/home/users/qd201969/gtopo30_0.9x1.25.nc")
phis=f0[2]
print(repr(phis))
phis=phis/9.8 # transfer from m2/s2 to m

# -------------- draw figure ---------------
titls=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    
cfp.setvars(file='month_'+da.attrs['long_name']+'.png')
cfp.gopen(figsize=[20, 20],rows=4,columns=3,wspace=0.1,hspace=0.015,bottom=0.5)
cfp.mapset(lonmin=lonl, lonmax=lonr, latmin=lats, latmax=latn)

for np in range(0,len(titls),1):#,len(f),1):
    #np = nm+1
    cfp.gpos(np+1)
    #cfp.levs(manual=[500,850])
    #cfp.cscale('/home/users/qd201969/uor_track/rwb.txt')
    cfp.levs(min=lev[0], max=lev[1], step=lev[2])
    cfp.cscale('precip2_17lev')
    cfp.con(var[np,:,:],x=ilon, y=ilat, ptype=1, 
            fill=True, lines=False, colorbar=None,title=da.attrs['long_name']+' '+titls[np])
    cfp.levs(manual=[1500,3000])
    cfp.con(phis,fill=False, lines=True, colors='k',linewidths=2)
    cfp.levs()

#cfp.levs(manual=[500,850])
#cfp.cscale('/home/users/qd201969/uor_track/rwb.txt')
cfp.cscale('precip2_17lev')
cfp.levs(min=lev[0], max=lev[1], step=lev[2])
cfp.cbar(position=[0.2, 0.48, 0.6, 0.01])
cfp.gclose()

subprocess.run('mogrify -bordercolor white -trim ./month_*.png',shell=True) 

