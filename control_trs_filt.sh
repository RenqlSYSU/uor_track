#!/bin/bash
cd ~/TRACK-1.5.2
pwd

lats=25
latn=45
lonl=60
lonr=80
option=0 #Genesis (0)/Lysis (1)/Passing(2)/Passing Time(3)/All Times(4)

echo "=========== trs filt (${lats}-${latn}N, ${lonl}-${lonr}E) ================"
file=~/ERA5-1HR/ff_trs_pos.oct-mar1979-2020_addmslp_addwind925_addwind10m
#file=~/ERA5-1HR/ff_trs_pos.apr-sep1979-2019_addmslp_addwind925_addwind10m
echo $file

utils/bin/box $file ${lats} ${latn} ${lonl} ${lonr} $option 0 0.0
mv ${file}.new ${file}.${option}_${lats}${latn}-${lonl}${lonr}
utils/bin/tr2nc ${file}.${option}_${lats}${latn}-${lonl}${lonr} s utils/TR2NC/tr2nc.meta

