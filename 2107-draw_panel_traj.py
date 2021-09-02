import cf
import cfplot as cfp
import numpy as np
import sys
import subprocess

f0=cf.read("/home/users/qd201969/gtopo30_0.9x1.25.nc")
phis=f0[2]
print(repr(phis))
phis=phis/9.8 # transfer from m2/s2 to m

path = '/home/users/qd201969/TRACK_result/'
filename = ['ff_trs_pos','tr_trs_pos.0_2550-4070']
lev = [850,500,300]
filt = False
nc = 1

cfp.setvars(file=filename[nc]+'.png')
cfp.gopen(rows=3, columns=1,hspace=0.25)
cfp.mapset(lonmin=-90, lonmax=150, latmin=10, latmax=70)

for nl in range(0,3,1):#,len(f),1):
    np=nl+1
    f=cf.read(path+'ERA5_VOR'+str(lev[nl])+'_6hr_2000_DET_T42/'+filename[nc]+'.nc')
    print(f)
    g = f[2]
    g = g*1e5
    print(g)

    if filt:
        x = g.aux('longitude').data
        y = g.aux('latitude').data
        tim = g.aux('time').data

        nt = 0
        ni = range(0,g.shape[0])
        nj = range(0,g.shape[1])

        for i in ni:
            inreg = 0
            h = g[i].where(g[i].array>1e+10, cf.masked)
            print(h.shape)
            h = np.squeeze(h.array)
            print(h.shape)

            for j in nj:
                if x[i,j] >= 30 and x[i,j] <= 110 and y[i,j] >= 20 and y[i,j] <= 55 and not x[i,j].mask and not y[i,j].mask:
                    jmax = np.argmax(h)
                    tmax = h[jmax]
                    longmax = x[i, jmax]
                    latmax  = y[i, jmax]
                    print(longmax, latmax, tmax)
                    inreg = 1
                    break

            if inreg == 1:
                if nt == 0:
                    tr_arb = g[i]
                    nt = 1
                else:
                    tr_arb = [tr_arb, g[i]]

        tr_arb = cf.aggregate(tr_arb, relaxed_identities=True)
        g = tr_arb[0]

    cfp.gpos(np)
    cfp.levs(manual=[1500,3000,4500])
    cfp.con(phis,fill=False, lines=True, colors='k',linewidths=2)
    cfp.levs()

    cfp.cscale('precip2_17lev')
    cfp.levs(min=0.0, max=8.0, step=0.5)
    cfp.traj(g, zorder=0, legend_lines=True, colorbar=False, linewidth=1.5, title=filename[nc]+' '+str(lev[nl]))

cfp.cbar(position=[0.9, 0.2, 0.01, 0.6], title='Relative Vorticity (Hz)*1e5',orientation='vertical')
cfp.gclose()
subprocess.run('mogrify -bordercolor white -trim /home/users/qd201969/ERA5-1HR/'+filename[nc]+'.png',shell=True) 

