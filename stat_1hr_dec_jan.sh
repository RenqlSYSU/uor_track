#!/bin/bash

cd ~/TRACK-1.5.2/
filname=$1
output=$2
season=$3 # 0 = monthly; 1 = seasonal
echo $filname
echo $output

if [ $season == 0 ];then
    nday=(31 28 31 30 31 30 31 31 30 31 30 31)
    #nmonth=(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec)
    nmonth=($(seq 1 1 12))
    frae=744
else
    nday=(90 92 92 91)
    nmonth=(DJF MAM JJA SON)
    frae=0
fi
echo ${nmonth[*]}
    
sed -i "21s/.*/${filname}/" indat/STATS.latlng_1hr.in
for nm in $(seq 1 1 ${#nday[*]});do #{1..${nm}};do
    fras=$((frae+1))
    frae=$((fras+24*nday[$((nm-1))]-1))
    echo "month=${nmonth[$((nm-1))]}, frame_s=${fras}, frame_e=${frae}"

    sed -i "34s/.*/${fras}/" indat/STATS.latlng_1hr.in
    sed -i "35s/.*/${frae}/" indat/STATS.latlng_1hr.in

    bin/track.linux < indat/STATS.latlng_1hr.in > ${output}_record
    mv outdat/stat_trs_scl.linux_1.nc ${output}_${nmonth[$((nm-1))]}.nc
done

