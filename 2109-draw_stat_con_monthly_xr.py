import cf
import cfplot as cfp
import subprocess
import matplotlib
import numpy.ma as ma
import sys
import xarray as xr
import numpy as np
import gc #garbage collector
import cartopy.crs as ccrs
matplotlib.use('Agg')

lonl=0  #0  #
lonr=150#360#
lats=15 #0  #
latn=70 #90 #

lev=[[0    ,320  ,20  ], # 0Feature Density
     [0    ,3.2  ,0.2 ], # 1Genesis Density
     [0    ,3.2  ,0.2 ], # 2Lysis Density
     [-1   ,1    ,1   ], # 3Mean Area
     [-0.8 ,0.8  ,0.1 ], # 4Mean Growth/Decay Rate
     [-1   ,1    ,1   ], # 5Mean Anisotropy
     [0    ,8    ,0.5 ], # 6Mean Lifetime
     [0    ,80   ,5   ], # 7Mean Speed
     [0    ,16   ,1   ], # 8Mean Intensity
     [-1.6 ,1.6  ,0.2 ], # 9Mean Tendency
     [-1   ,1    ,1   ], # 10Spare1
     [-1   ,1    ,1   ], # 11Spare2
     [0    ,80   ,5   ], # 12Std of Speed
     [0    ,3.2  ,0.2 ], # 13Std of Intensity
     [0    ,16   ,1   ], # 14Track Density
     [-1   ,1    ,1   ], # 15X-component of Mean Orientation Vector
     [-40  ,40   ,5   ], # 16X-component of Mean Velocity
     [-1   ,1    ,1   ], # 17Y-component of Mean Orientation Vector
     [-40  ,40   ,5  ]] # 18Y-component of Mean Velocity

draw_var = ["fden","gden","lden","marea","mgdr","",
            "mlif","msp" ,"mstr","mten" ,""    ,"",
            ""    ,""    ,"tden",""    ] # 7 variables
#draw=[8,9,6]
#draw=[1,2,14]
draw=[14]
#draw=[1,2,14,8,9,6]

if int(sys.argv[1]) == 0:
    filename = 'ff_250_1980-2020_2_2545-6080'
    files = '/home/users/qd201969/ERA5-1HR-lev/statistic/'+filename+'_stat.nc'
    level = 2
    figtitle = '250'
    dbox = 0 
else:
    filename = sys.argv[1] #'ff_250_500_no'
    files = sys.argv[2] #'/home/users/qd201969/ERA5-1HR-lev/match'+filt+'/statistic/'+filename+'_stat_'
    level = int(sys.argv[3])
    figtitle = sys.argv[4]
    dbox = int(sys.argv[5])

if dbox >= 1 :
    flats = int(sys.argv[6])
    flatn = int(sys.argv[7])
    flonl = int(sys.argv[8])
    flonr = int(sys.argv[9])

if level == 0: # use total level bar 
    lev[ 0][:]=[0    ,1600 ,100 ]
    lev[ 1][:]=[0    ,8    ,0.5 ]
    lev[ 2][:]=[0    ,8    ,0.5 ]
    lev[14][:]=[0    ,48   ,3   ]
if level == 2: # 
    lev[ 1][:]=[0    ,4.8  ,0.3 ]
    lev[ 2][:]=[0    ,4.8  ,0.3 ]
    lev[14][:]=[0    ,24   ,1.5 ]
titls=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

f = xr.open_dataset(files)
lat = f.lat
lon = f.long
ilon = lon[(lon>=lonl) & (lon<=lonr)]
ilat = lat[(lat>=lats) & (lat<=latn)]

ds = xr.open_dataset('/home/users/qd201969/data/ERA5_mon_u_1979-2020.nc')
da = ds['u'].sel(level=200,expver=1,longitude=ilon,latitude=ilat,method="nearest").load()
# increased performance by loading data into memory first, e.g., with load()
uwnd = da.groupby(da.time.dt.month).mean('time')
print(uwnd)
del ds, da
gc.collect()

ds = xr.open_dataset("/home/users/qd201969/gtopo30_0.9x1.25.nc")
phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").load()
print(repr(phis))
phis = phis/9.8 # transfer from m2/s2 to m
del ds
gc.collect()

for nv in range(0,len(draw),1):#,len(f),1):
    var = f[draw_var[draw[nv]]].sel(long=ilon,lat=ilat).load()
    if draw[nv] == 9:
        var=var*24
    
    if draw[nv] > 2 and draw[nv] != 14:
        tden = f['tden'].sel(long=ilon,lat=ilat).load()
        mask = tden < 1.0
        var.values=ma.array(var.values,mask=mask)
    
    cfp.setvars(file='month_'+filename+var.long_name+'.png')
    cfp.gopen(figsize=[20, 20],rows=4,columns=3,wspace=0.1,hspace=0.015,bottom=0.5)
    cfp.mapset(lonmin=lonl, lonmax=lonr, latmin=lats, latmax=latn)
    
    for nm in range(0,len(titls),1):#,len(f),1):
        np = nm+1
        cfp.gpos(np)
        cfp.levs(min=lev[draw[nv]][0], max=lev[draw[nv]][1], step=lev[draw[nv]][2])
        if lev[draw[nv]][0] < 0 :
            cfp.cscale('BlueDarkRed18')
        else:
            cfp.cscale('precip2_17lev')
        cfp.con(var[nm,:,:],x=ilon, y=ilat, ptype=1,fill=True,
                lines=False,colorbar=None,title=figtitle+' '+titls[nm])
        cfp.levs(manual=[1500,3000]) # topography 
        cfp.con(phis,x=ilon, y=ilat, ptype=1,fill=False,
                lines=True, colors='k',linewidths=2)
        cfp.levs(manual=[30,32]) # jet stream
        cfp.con(uwnd[nm,:,:],x=ilon, y=ilat, ptype=1, fill=False,
                lines=True, colors='blueviolet',linewidths=2)
        if dbox >= 1 :
            cfp.plotvars.mymap.plot([flonl,flonl,flonr,flonr,flonl],[flatn,flats,flats,flatn,flatn], 
                    linewidth=4, color='grey', transform=ccrs.PlateCarree()) # filter box

    if lev[draw[nv]][0] < 0 :
        cfp.cscale('BlueDarkRed18')
    else:
        cfp.cscale('precip2_17lev')
    cfp.levs(min=lev[draw[nv]][0], max=lev[draw[nv]][1], step=lev[draw[nv]][2])
    cfp.cbar(position=[0.2, 0.48, 0.6, 0.01], title=var.long_name)
    cfp.gclose()

#subprocess.run('mogrify -bordercolor white -trim ./month_*.png',shell=True) 
#subprocess.run('mogrify -bordercolor white -trim /home/users/qd201969/ERA5-1HR/month_*.png',shell=True) 

