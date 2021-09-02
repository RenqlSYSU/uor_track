#!/bin/bash

nmonth=$(seq 4 1 9)
nday=(31 28 31 30 31 30 31 31 30 31 30 31)

fras=0
frae=0
for nm in ${nmonth[*]};do
    fras=$((frae+1))
    frae=$((fras+24*nday[$((nm-1))]-1))
    echo "month=${nm},frame_s=${fras},frame_e=${frae}"
done

for ny in {1983..2020};do
    echo $ny
done

