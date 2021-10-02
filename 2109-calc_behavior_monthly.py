#!/usr/bin/env python
'''
read total cyclone number in ff_250_1980-2020_2_3045-5960
then use box to filter different behaviors cyclones
plot the filter box in the trajectory figure

20210928
'''

import xarray as xr
import numpy as np
import sys, os, subprocess, linecache
import matplotlib.pyplot as plt

def calc_month( frae, month, nday, filname):
    var  = np.empty( [len(month)], dtype=int ) # 4 or 12 
    for nm in range(0,len(month),1):
        fras = frae+1
        frae = fras+24*nday[nm]-1
        print("month=%s, frame_s=%d, frame_e=%d" %(month[nm],fras,frae))
        rf = open("/home/users/qd201969/TRACK-1.5.2/region.dat","w")
        rf.write("0 360\n-90 90\n"+str(fras)+" "+str(frae))
        rf.close()

        ret=subprocess.Popen("utils/bin/censemble2 "+filname+" "\
                +filname+" 0 100 10 1 0 0 0 0 s 0 1 > rec",shell=True)
        ret.wait()

        a=linecache.getline("/home/users/qd201969/TRACK-1.5.2/match_ens1_yes.dat", 4).strip().split(" ",1)
        term = a[1].strip().split(" ",1)
        var[nm] = int(term[0])
        linecache.clearcache()
    return var

if sys.argv[1] == "-1" :
    option=2 #int(sys.argv[1]) #Genesis (0)/Lysis (1)/Passing(2)/Passing Time(3)/All Times(4)
    flats = 30 #int(sys.argv[2])
    flatn = 45 #int(sys.argv[3])
    flonl = 59.9 #int(sys.argv[4])
    flonr = 60 #int(sys.argv[5])
    prefix = "ff"
    season = 0 # 0 monthly, 1 seasonal
else:
    option= int(sys.argv[1]) 
    flats = int(sys.argv[2])
    flatn = int(sys.argv[3])
    flonl = int(sys.argv[4])
    flonr = int(sys.argv[5])
    prefix = int(sys.argv[6])
    season = int(sys.argv[7])
suffix=str(option)+"_"+str(flats)+str(flatn)+"-"+str(flonl)+str(flonr)
figdir = "/home/users/qd201969/uor_track/fig/behv_month_"+suffix
fileout="/home/users/qd201969/uor_track/mdata/behv_month_"+suffix+".nc"
calcbehv = 0

behv = ["ALL"  ,"NTN"  ,"STN"  ,"PAS"  ,"LYS"  ,"DIF"]
lats = [flats  ,flatn  ,flats-1,flats  ,flats  ]
latn = [flatn  ,flatn+1,flats  ,flatn  ,flatn  ]
lonl = [flonl  ,flonl  ,flonl  ,90     ,flonr  ]
lonr = [flonr  ,105    ,105    ,105    ,90     ]
opti = [option ,2      ,2      ,2      ,1      ]
lev = [850,500,250]

if season == 0:
    nday=[365,31,28,31,30,31,30,31,31,30,31,30,31]
    month=["ALL","Jan","Feb","Mar","Apr","May","Jun",\
            "Jul","Aug","Sep","Oct","Nov","Dec"]
    frae=744
else:
    nday=[365,90,92,92,91]
    month=["ALL","DJF","MAM","JJA","SON"]
    frae=0

imonth=range(1,len(month),1)
var  = np.empty( [len(month),len(lev),len(behv)],dtype=int )  
perc = np.empty( [len(month),len(lev),len(behv)],dtype=float )  
path = '/home/users/qd201969/ERA5-1HR-lev/'
os.chdir("/home/users/qd201969/TRACK-1.5.2")
#===============================================
# box filer cyclone, read cyclone number
#===================================================
if calcbehv == 1:
    for nl in range(0,len(lev),1):
        filname  = path+prefix+"_"+str(lev[nl])+"_1980-2020_"+suffix

        if not os.path.isfile(filname) :
            filname0 = path+prefix+"_"+str(lev[nl])+"_1980-2020"
            ret=subprocess.Popen("utils/bin/box "+filname0+" "+str(lats[0])+" "+str(latn[0])+\
                    " "+str(lonl[0])+" "+str(lonr[0])+" "+str(opti[0])+" 0 0.0",shell=True)
            ret.wait()
            subprocess.run("mv "+filname0+".new "+filname,shell=True)

        a=linecache.getline(filname, 4).strip() #读取特定行数据
        a=a.split(" ",1)
        term = a[1].strip().split(" ",1)
        var[0,nl,0] = int(term[0])
        linecache.clearcache()
        var[1:len(month),nl,0] = calc_month( frae, month[1:len(month)], nday[1:len(month)], filname)

        for nr in range(1,len(lats),1):
            ret=subprocess.Popen("utils/bin/box "+filname+" "+str(lats[nr])+" "+str(latn[nr])+\
                    " "+str(lonl[nr])+" "+str(lonr[nr])+" "+str(opti[nr])+" 0 0.0",shell=True)
            ret.wait()
            a=linecache.getline(filname+".new", 4).strip().split(" ",1)
            term = a[1].strip().split(" ",1)
            var[0,nl,nr] = int(term[0])
            linecache.clearcache()
            var[1:len(month),nl,nr] = calc_month( frae, month[1:len(month)], nday[1:len(month)], filname+".new")
            subprocess.run("rm "+filname+".new",shell=True)

    var[:,:,-1] = var[:,:,0]-np.sum(var[:,:,1:5],axis=2)

    f = open("/home/users/qd201969/uor_track/behv-monthly-"+prefix+"-"+suffix,'w')
    f.write(path+prefix+"_1980-2020_"+suffix+"\n")
    f.write("*"*50+"\n")
    f.write("lats"+str(lats)+"\n")
    f.write("latn"+str(latn)+"\n")
    f.write("lonl"+str(lonl)+"\n")
    f.write("lonr"+str(lonr)+"\n")
    f.write("opti"+str(opti)+"\n")
    f.write("*"*50+"\n")
    for nm in range(0,len(month),1):
        f.write("\n")
        f.write("* "+month[nm]+" number (percent)\n")
        f.write("* lev, "+str(behv).strip('[').strip(']').replace('\'','')+"\n")
        for nl in range(0,len(lev),1):
            perc[nm,nl,:] = 100*var[nm,nl,:]/var[nm,nl,0]
            output = str(var[nm,nl,0])
            for nr in range(1,len(lats)+1,1):
                output = output+", "+str(var[nm,nl,nr])+" ("+str(np.around(perc[nm,nl,nr],decimals=2))+")"
            f.write("* "+str(lev[nl])+", "+output+"\n")
    f.close()

    ds = xr.Dataset(
            {
                "num" : (["month", "lev", "behv"], var),
                "perc": (["month", "lev", "behv"], perc),
                },
            coords={
                "month": range(0,len(month),1), 
                "lev"  : (["lev"], lev),
                "behv" : (["behv"],behv),
                },
            )
    ds.attrs["description"] = "number and percent of cycle"
    ds.to_netcdf(fileout,"w")

#===============================================
# draw figure 
#===================================================
nrow = 3
ncol = 2
bmlo = 0.4
BIGFONT=22
MIDFONT=14
SMFONT=10
if calcbehv != 1:
    fvar = xr.open_dataset(fileout)
    var = fvar['num'].data
    perc= fvar['perc'].data

fig = plt.figure(figsize=(9,9),dpi=200)
ax = fig.subplots(nrow, ncol) #sharex=True, sharey=True
for nl in range(0,len(lev),1):
    ax[nl][0].set_title(str(lev[nl])+" number "+suffix,fontsize=MIDFONT)
    ax[nl][1].set_title(str(lev[nl])+" percent "+suffix,fontsize=MIDFONT)
    for nr in range(1,len(lats),1):
        ax[nl][0].plot(imonth,var[1:len(month),nl,nr])
        ax[nl][1].plot(imonth,perc[1:len(month),nl,nr])
    ax[nl][0].legend(behv[1:5], loc='upper left')
    ax[nl][1].legend(behv[1:5], loc='upper left')

fig.tight_layout(w_pad=0.5,h_pad=1) #,rect=(0,bmlo,1,1)
fig.savefig(figdir+".png")
#fig.savefig(figdir+".png", bbox_inches='tight',pad_inches=0.01)

print(path+prefix+"_1980-2020_"+suffix)
print("lats"+str(lats))
print("latn"+str(latn))
print("lonl"+str(lonl))
print("lonr"+str(lonr))
print("opti"+str(opti))
print("*"*50)
print("lev "+str(behv).strip('[').strip(']').replace(',','').replace('\'',''))
for nl in range(0,len(lev),1):
    print(str(lev[nl])+" "+str(var[0,nl,:]).strip('[').strip(']'))

print("")
print("lev "+str(behv).strip('[').strip(']').replace(',','').replace('\'',''))
for nl in range(0,len(lev),1):
    print(str(lev[nl])+" "+str(perc[0,nl,:]).strip('[').strip(']'))
