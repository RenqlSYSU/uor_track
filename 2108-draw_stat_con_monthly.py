import cf
import cfplot as cfp
import subprocess
import matplotlib
import numpy.ma as ma
import sys
matplotlib.use('Agg')

lonl=0  #0  #
lonr=360#150#
lats=0  #15 #
latn=90 #70 #

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

draw=[8,9,6]
#draw=[1,2,14]
#draw=[1,2,14,8,9,6]
level = int(sys.argv[3])
if level == 0: # use total level bar 
    lev[ 0][:]=[0    ,1600 ,100 ]
    lev[ 1][:]=[0    ,8    ,0.5 ]
    lev[ 2][:]=[0    ,8    ,0.5 ]
    lev[14][:]=[0    ,48   ,3   ]

if level == 2: # 
    lev[ 1][:]=[0    ,4.8  ,0.3 ]
    lev[ 2][:]=[0    ,4.8  ,0.3 ]
    lev[14][:]=[0    ,24   ,1.5 ]

f0=cf.read("/home/users/qd201969/gtopo30_0.9x1.25.nc")
phis=f0[2]
print(repr(phis))
phis=phis/9.8 # transfer from m2/s2 to m

filename = sys.argv[1] #'ff_250_500_no'
files = sys.argv[2] #'/home/users/qd201969/ERA5-1HR-lev/match'+filt+'/statistic/'+filename+'_stat_'
figname = sys.argv[4]
#files='/home/users/qd201969/ERA5-1HR-lev/statistic/'+filename+'_stat_'
#files='/home/users/qd201969/ERA5-1HR/stat_clim_'+filt
titls=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

f1   = cf.read(files+'1.nc')
for nv in range(0,len(draw),1):#,len(f),1):
    g1=f1[draw[nv]]
    
    cfp.setvars(file='month_'+filename+g1.long_name+'.png')
    cfp.gopen(figsize=[20, 20],rows=4,columns=3,wspace=0.1,hspace=0.015,bottom=0.5)
    cfp.mapset(lonmin=lonl, lonmax=lonr, latmin=lats, latmax=latn)
    
    for nm in range(0,len(titls),1):#,len(f),1):
        np = nm+1
        f=cf.read(files+str(np)+'.nc')
        g=f[draw[nv]] # read data
        if draw[nv] == 9:
            g=g*24
        
        if draw[nv] > 2 and draw[nv] != 14:
            trd  = f[14]
            mask = trd < 1.0
            g.data=ma.array(g,mask=mask)

        cfp.gpos(np)
        cfp.levs(min=lev[draw[nv]][0], max=lev[draw[nv]][1], step=lev[draw[nv]][2])
        if lev[draw[nv]][0] < 0 :
            cfp.cscale('BlueDarkRed18')
        else:
            cfp.cscale('precip2_17lev')

        cfp.con(g,fill=True,lines=False,colorbar=None,title=figname+' '+titls[nm])
        cfp.levs(manual=[1500,3000])
        cfp.con(phis,fill=False, lines=True, colors='k',linewidths=2)
        cfp.levs()
    
    if lev[draw[nv]][0] < 0 :
        cfp.cscale('BlueDarkRed18')
    else:
        cfp.cscale('precip2_17lev')
    cfp.levs(min=lev[draw[nv]][0], max=lev[draw[nv]][1], step=lev[draw[nv]][2])
    cfp.cbar(position=[0.2, 0.48, 0.6, 0.01], title=g1.long_name)
    cfp.gclose()

#subprocess.run('mogrify -bordercolor white -trim ./month_*.png',shell=True) 
#subprocess.run('mogrify -bordercolor white -trim /home/users/qd201969/ERA5-1HR/month_*.png',shell=True) 

