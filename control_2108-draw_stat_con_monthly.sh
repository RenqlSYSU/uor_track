#!/bin/bash

lev=(850 500 250)
OUTDIR=/home/users/qd201969/ERA5-1HR-lev/
prefix=ff # tr is original cyclone ; ff is the filtered cyclone

filt=1 #filt=1, then use filter track to match

# filter area
lats=25
latn=45
lonl=60
lonr=80
option=2 #Genesis (0)/Lysis (1)/Passing(2)/Passing Time(3)/All Times(4)
if [ $filt == 1 ]; then 
    suffix=_${option}_${lats}${latn}-${lonl}${lonr}
else
    suffix=""
fi

figname=('only_250' '250-500' '250-500-850' \
         'only_500' '500-850' '500-250' '500-250-850' \
         'only_850' '850-500' '850-500-250')
cd ${OUTDIR}match${suffix}/
level=0 # 1=total cyclone level, 0=filt cyclone level 
path=$(pwd)

np=0
for filename in ${prefix}_*;do
    echo ${filename}
    file=${path}/statistic/${filename}_stat_
    python ~/uor_track/2108-draw_stat_con_monthly.py \
        ${suffix} ${filename} ${file} ${level} ${figname[$np]}
    np=$((np+1))
done

mogrify -bordercolor white -trim ./month_*.png

