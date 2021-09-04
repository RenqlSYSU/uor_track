#!/bin/bash

nmonth=$(seq 1 1 12)
nday=(31 28 31 30 31 30 31 31 30 31 30 31)

fras=0
frae=0
for nm in ${nmonth[*]};do
    fras=$((frae+1))
    frae=$((fras+24*nday[$((nm-1))]-1))
    echo "month=${nm},frame_s=${fras},frame_e=${frae}"
done

lev=(850 500 250)
OUTDIR=/home/users/qd201969/ERA5-1HR-lev/
prefix=ff # tr is original cyclone ; ff is the filtered cyclone

# filter area
lats=25
latn=45
lonl=60
lonr=80
option=2 #Genesis (0)/Lysis (1)/Passing(2)/Passing Time(3)/All Times(4)
suffix=_${option}_${lats}${latn}-${lonl}${lonr}

cd ${OUTDIR}match${suffix}
for file in ${prefix}_* ; do
    cd ~/TRACK-1.5.2
    echo $file
done

