#!/usr/bin/env python
'''
read total cyclone number in ff_250_1980-2020_2_3045-5960
then use box to filter different behaviors cyclones
plot the filter box in the trajectory figure

20210928
'''

import cf
import cfplot as cfp
import numpy as np
import sys, os, subprocess, linecache
import cartopy.crs as ccrs

if len(sys.argv) < 2 :
    option=2 #int(sys.argv[1]) #Genesis (0)/Lysis (1)/Passing(2)/Passing Time(3)/All Times(4)
    flats = 0 #int(sys.argv[2])
    flatn = 90 #int(sys.argv[3])
    flonl = 60 #int(sys.argv[4])
    flonr = 60 #int(sys.argv[5])
    prefix = "ff"
else:
    option= int(sys.argv[1]) 
    flats = int(sys.argv[2])
    flatn = int(sys.argv[3])
    flonl = int(sys.argv[4])
    flonr = int(sys.argv[5])
    prefix = int(sys.argv[6])
suffix=str(option)+"_"+str(flats)+str(flatn)+"-"+str(flonl)+str(flonr)

flonr2 = 90
behv = ["ALL" ,"NTN" ,"STN" ,"PAS" ,"LYS" ]#,"DIF"]
lats = [flats ,flatn ,flats ,flats ,flats ]
latn = [flatn ,flatn ,flats ,flatn ,flatn ]
lonl = [flonl ,flonl ,flonl ,flonr2,flonr ]
lonr = [flonr ,flonr2,flonr2,flonr2,flonr2]
opti = [option ,2      ,2      ,2      ,1      ]
lev = [850,500,250]
var  = np.empty( [len(lev),len(behv)],dtype=int )  
perc = np.empty( [len(lev),len(behv)],dtype=float )  

path = '/home/users/qd201969/ERA5-1HR-lev/'
drawbox = 1
calcbehv = 0

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
        nm = a[1].strip().split(" ",1)
        var[nl,0] = int(nm[0])
        linecache.clearcache()
        
        for nr in range(1,len(lats),1):
            ret=subprocess.Popen("utils/bin/box "+filname+" "+str(lats[nr])+" "+str(latn[nr])+\
                    " "+str(lonl[nr])+" "+str(lonr[nr])+" "+str(opti[nr])+" 0 0.0",shell=True)
            ret.wait()
            a=linecache.getline(filname+".new", 4).strip().split(" ",1)
            nm = a[1].strip().split(" ",1)
            var[nl,nr] = int(nm[0])
            linecache.clearcache()
            subprocess.run("rm "+filname+".new",shell=True)

        var[nl,-1] = var[nl,0]-np.sum(var[nl,1:5])
        perc[nl,:] = 100*var[nl,:]/var[nl,0]

    f = open("/home/users/qd201969/uor_track/behv-"+prefix,'a')
    f.write("\n\n")
    f.write(path+prefix+"_1980-2020_"+suffix+"\n")
    f.write("lats"+str(lats)+"\n")
    f.write("latn"+str(latn)+"\n")
    f.write("lonl"+str(lonl)+"\n")
    f.write("lonr"+str(lonr)+"\n")
    f.write("opti"+str(opti)+"\n")
    f.write("*"*50+"\n")
    f.write("lev "+str(behv).strip('[').strip(']').replace(',','').replace('\'','')+"\n")
    for nl in range(0,len(lev),1):
        f.writelines(str(lev[nl])+" "+str(var[nl,:]).strip('[').strip(']')+"\n")

    f.write("\n")
    f.write("lev "+str(behv).strip('[').strip(']').replace(',','').replace('\'','')+"\n")
    for nl in range(0,len(lev),1):
        f.write(str(lev[nl])+" "+str(np.around(perc[nl,:],decimals=2)).strip('[').strip(']')+"\n")
    f.close()
    
    print(path+prefix+"_1980-2020_"+suffix)
    print("lats"+str(lats))
    print("latn"+str(latn))
    print("lonl"+str(lonl))
    print("lonr"+str(lonr))
    print("opti"+str(opti))
    print("*"*50)
    print("lev "+str(behv).strip('[').strip(']').replace(',','').replace('\'',''))
    for nl in range(0,len(lev),1):
        print(str(lev[nl])+" "+str(var[nl,:]).strip('[').strip(']'))

    print("")
    print("lev "+str(behv).strip('[').strip(']').replace(',','').replace('\'',''))
    for nl in range(0,len(lev),1):
        print(str(lev[nl])+" "+str(perc[nl,:]).strip('[').strip(']'))

#===============================================
# draw figure 
#===================================================
if drawbox == 1:
    os.chdir("/home/users/qd201969/uor_track/fig")
    f0=cf.read("/home/users/qd201969/gtopo30_0.9x1.25.nc")
    phis=f0[2]
    print(repr(phis))
    phis=phis/9.8 # transfer from m2/s2 to m
    
    for nr in range(0,len(lats)-1,1):
        if nr == 0:
            suffix2 = ""
        else:
            suffix2 = "_"+str(opti[nr])+"_"+str(lats[nr])+str(latn[nr])+"-"+str(lonl[nr])+str(lonr[nr])
        
        cfp.setvars(file="traj-"+prefix+"_"+suffix+suffix2+".png")
        cfp.gopen(rows=3, columns=1 ,hspace=0.25)#,bottom=0.2
        cfp.mapset(lonmin=0, lonmax=150, latmin=10, latmax=70)
        for nl in range(0,3,1):#,len(f),1):
            np=nl+1
            
            cfp.gpos(np)
            cfp.levs(manual=[1500,3000,4500])
            cfp.con(phis,fill=False, lines=True, colors='k',linewidths=2,\
                    title=prefix+' '+suffix+" "+str(lev[nl])+" "+behv[nr])
            cfp.levs()
            
            for nrr in range(0,len(lats),1):
                cfp.plotvars.mymap.plot([lonl[nrr],lonl[nrr],lonr[nrr],lonr[nrr],lonl[nrr]],[latn[nrr],lats[nrr],lats[nrr],latn[nrr],latn[nrr]], 
                        linewidth=4, color='k', transform=ccrs.PlateCarree()) # filter box

            if behv[nr] == "PAS" and lev[nl] == 850:
                continue
            
            filname  = path+prefix+"_"+str(lev[nl])+"_1980-2020_"+suffix+suffix2
            if not os.path.isfile(filname+'.nc') :
                ret=subprocess.Popen("/home/users/qd201969/TRACK-1.5.2/utils/bin/tr2nc \
                        "+filname+" s /home/users/qd201969/TRACK-1.5.2/utils/TR2NC/tr2nc.meta",shell=True)
                ret.wait()
            f=cf.read(filname+'.nc')
            print(f)
            g = f[2]
            g = g*1e5
            print(g)

            cfp.cscale('precip2_17lev')
            cfp.levs(min=0.0, max=8.0, step=0.5)
            cfp.traj(g, zorder=0, legend_lines=True, colorbar=False, linewidth=1.5)
            
        cfp.cbar(position=[0.75, 0.2, 0.01, 0.6], title='Relative Vorticity (Hz)*1e5',orientation='vertical')
        cfp.gclose()
    subprocess.run("mogrify -bordercolor white -trim ./traj-"+prefix+"_"+suffix+"*.png",shell=True) 

