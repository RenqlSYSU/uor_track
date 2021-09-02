#!/bin/bash
cd ~/TRACK-1.5.2
pwd

lats=25
latn=50
lonl=40
lonr=70
option=0 #Genesis (0)/Lysis (1)/Passing(2)/Passing Time(3)/All Times(4)
filt=1 # 1 mean to do the region filter

echo "=========== trs filt (${lats}-${latn}N, ${lonl}-${lonr}E) ================"
file[0]=~/ERA5-1HR/ff_trs_pos.oct-mar1979-2020_addmslp_addwind925_addwind10m
file[1]=~/ERA5-1HR/ff_trs_pos.apr-sep1979-2019_addmslp_addwind925_addwind10m
outfile=~/ERA5-1HR/stat_clim
filname[0]="\/home\/users\/qd201969\/ERA5-1HR\/ff_trs_pos.oct-mar1979-2020_addmslp_addwind925_addwind10m"
filname[1]="\/home\/users\/qd201969\/ERA5-1HR\/ff_trs_pos.apr-sep1979-2019_addmslp_addwind925_addwind10m"

nday=(31 29 31 30 31 30 31 31 30 31 30 31)
echo ${#nday[*]} #print the length of nday

if [ ${filt} == 1 ]; then
    outfile=${outfile}_${option}_${lats}${latn}-${lonl}${lonr}
fi

for nc in {0..1};do
    echo '' 
    echo ${file[$nc]}
    if [ ${filt} == 1 ]; then
        utils/bin/box ${file[$nc]} ${lats} ${latn} ${lonl} ${lonr} $option 0 0.0 << EOF
n
EOF
        mv ${file[$nc]}.new ${file[$nc]}.${option}_${lats}${latn}-${lonl}${lonr}
        utils/bin/tr2nc ${file[$nc]}.${option}_${lats}${latn}-${lonl}${lonr} s ~/ERA5-1HR/tr2nc.meta 
        term=${filname[$nc]}.${option}_${lats}${latn}-${lonl}${lonr}
    fi

    echo ${term}
    sed -i "21s/.*/${term}/" indat/STATS.latlng_addf_1hr.in
    sed -i "34s/.*/1/" indat/STATS.latlng_addf_1hr.in
    sed -i "35s/.*/100000000/" indat/STATS.latlng_addf_1hr.in
    #bin/track.linux < indat/STATS.latlng_addf_1hr.in >> record4
    
    if [ ${nc} -eq 0 ]
    then 
        nmonth=(10 11 12 1 2 3)
    else
        nmonth=$(seq 4 1 9)
        #echo "Apr-Sep STATS================"
        #mv outdat/stat_trs_scl.linux_1.nc ${outfile}_apr-sep1979-2019.nc
    fi
    fras=0
    frae=0
    
    echo ${nmonth[*]}
    for nm in ${nmonth[*]};do
        fras=$((frae+1))
        frae=$((fras+24*nday[$((nm-1))]-1))
        echo "month=${nm}, frame_s=${fras}, frame_e=${frae}"

        sed -i "34s/.*/${fras}/" indat/STATS.latlng_addf_1hr.in
        sed -i "35s/.*/${frae}/" indat/STATS.latlng_addf_1hr.in

        bin/track.linux < indat/STATS.latlng_addf_1hr.in >> record4
        mv outdat/stat_trs_scl.linux_1.nc ${outfile}_${nm}.nc
    done

    unset nmonth
done

