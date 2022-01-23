#!/bin/bash

cd /gws/nopw/j04/ncas_generic/users/renql/ERA5_hourly/z/

lev=(250 500 850)

#for ny in {1979..2006};do
#    if [ ! -f ERA5_250z_${ny}.nc ];then
#        for nl in ${lev[@]};do
#            cdo sellevel,${nl} ERA5_NH_z_${ny}.nc ERA5_${nl}z_${ny}.nc
#        done
#        rm ERA5_NH_z_${ny}.nc
#    fi
#done

for ny in {1982..2006};do
ny1=$((ny+1))
ny0=$((ny-1))
for nl in ${lev[@]};do
    if [ ! -f ERA5_${nl}z_${ny1}-01.nc ];then
        cdo selmonth,1 ERA5_${nl}z_${ny1}.nc ERA5_${nl}z_${ny1}-01.nc
        cdo selmonth,12 ERA5_${nl}z_${ny0}.nc ERA5_${nl}z_${ny0}-12.nc
    fi
    
    if [ ! -f ERA5_z${nl}_1hr_dec-jan${ny}.nc ];then
        cdo -b F32 -mergetime ERA5_${nl}z_${ny0}-12.nc ERA5_${nl}z_${ny}.nc \
            ERA5_${nl}z_${ny1}-01.nc ERA5_z${nl}_1hr_dec-jan${ny}.nc   
    fi 
    rm ERA5_${nl}z_${ny1}-01.nc
    rm ERA5_${nl}z_${ny0}-12.nc
    rm ERA5_${nl}z_${ny0}.nc
done
done

