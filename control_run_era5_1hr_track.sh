#!/bin/bash
#SBATCH --partition=cluster
#SBATCH --output=/home/users/qd201969/ERA5-1HR-lev/TMP/slurm-%j.out 
#SBATCH --error=/home/users/qd201969/ERA5-1HR-lev/TMP/slurm-%j.err
#SBATCH --job-name=ew
#SBATCH --time=10-00:00:00

lev=(850 500 250)
OUTDIR=/home/users/qd201969/ERA5-1HR-lev/
file=ff_trs_pos

pro=4  # 0 = run track ; 1 = combine ; 2 = statistics ; 3 = trs filt ; 4 = match
filt=0 #if pro=4, filt=1, then use filter track to match
nl1=2
nl2=1

# filter area
lats=25
latn=50
lonl=40
lonr=70
option=0 #Genesis (0)/Lysis (1)/Passing(2)/Passing Time(3)/All Times(4)

if [ $pro == 0 ];then
    nl=2
    cd ~/TRACK-1.5.2/
    pwd
    echo "run track run_era5_1hr_qiaoling.csh ${lev[$nl]}"
    for ny in {2007..2020};do
        ./run_era5_1hr_qiaoling.csh ${ny} ${lev[$nl]} > ~/ERA5-1HR-lev/record_${lev[$nl]}_${ny}
    done
fi

if [ $pro == 1 ];then
    cd $OUTDIR
    pwd
    echo "combine ${file}"
    for nl in {1..1};do
        echo 41 > combine.in_${lev[$nl]}
        echo 1 >> combine.in_${lev[$nl]}
        
        for ny in {1980..2020};do
            echo ${OUTDIR}ERA5_VOR${lev[$nl]}_1hr_${ny}_DET/${file} >> combine.in_${lev[$nl]}
        
            if [ ! -f "${OUTDIR}ERA5_VOR${lev[$nl]}_1hr_${ny}_DET/${file}" ]; then
                echo "${OUTDIR}ERA5_VOR${lev[$nl]}_1hr_${ny}_DET/${file} does not exist"
                #rm -irf ${OUTDIR}ERA5_VOR${lev[$nl]}_1hr_${ny}_DET
                #~/TRACK-1.5.2/run_era5_1hr_qiaoling.csh ${ny} ${lev[$nl]}
            else
                #awk '{if(NF==4 && ($2 > 400 || $3 > 100 || $4 > 1000)) print FILENAME, $0}' \
                #    ${OUTDIR}ERA5_VOR${lev[$nl]}_1hr_${ny}_DET/${file}
                num=$(awk '{if(NF==4 && ($2 > 400 || $3 > 100 || $4 > 1000)) print FILENAME, $0}' ${OUTDIR}ERA5_VOR${lev[$nl]}_1hr_${ny}_DET/${file} | awk 'END{print NR}')
                if [ ${num} -ge 2 ]; then 
                    echo "Need to rerun ${OUTDIR}ERA5_VOR${lev[$nl]}_1hr_${ny}_DET, whose error line is ${num}"
                    #rm -irf ${OUTDIR}ERA5_VOR${lev[$nl]}_1hr_${ny}_DET
                    #~/TRACK-1.5.2/run_era5_1hr_qiaoling.csh ${ny} ${lev[$nl]}
                fi
            fi
        done
        
        ~/TRACK-1.5.2/utils/bin/combine < combine.in_${lev[$nl]}
        mv ./combined_tr_trs ./${file}_${lev[$nl]}_1980-2020 
    done
fi

if [ $pro == 2 ];then
    echo "=========== statistics ================"
    cd ~/TRACK-1.5.2/
    pwd
    for nl in {2..2};do
        file=ff_trs_pos_${lev[$nl]}_1980-2020
        filname="\/home\/users\/qd201969\/ERA5-1HR-lev\/${file}"
        echo $filname

        for season in {0..1};do # 0 = monthly; 1 = seasonal
        if [ $season == 0 ];then
            nday=(31 28 31 30 31 30 31 31 30 31 30 31)
            nmonth=(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec)
            frae=744
        else
            nday=(90 92 92 91)
            nmonth=(DJF MAM JJA SON)
            frae=0
        fi
            
        sed -i "21s/.*/${filname}/" indat/STATS.latlng_1hr.in
        for nm in $(seq 1 1 ${#nday[*]});do #{1..${nm}};do
            fras=$((frae+1))
            frae=$((fras+24*nday[$((nm-1))]-1))
            echo "month=${nmonth[$((nm-1))]}, frame_s=${fras}, frame_e=${frae}"

            sed -i "34s/.*/${fras}/" indat/STATS.latlng_1hr.in
            sed -i "35s/.*/${frae}/" indat/STATS.latlng_1hr.in

            bin/track.linux < indat/STATS.latlng_1hr.in > record4_${file}
            mv outdat/stat_trs_scl.linux_1.nc ${OUTDIR}${file}_stat_${nmonth[$((nm-1))]}.nc
        done
        done
    done
fi

if [ $pro == 3 ];then
    echo "=========== trs filt (${lats}-${latn}N, ${lonl}-${lonr}E) ================"
    cd ~/TRACK-1.5.2/
    pwd
    for nl in {0..2};do
        echo "lev: ${lev[$nl]}"
        file=~/TRACK_result/ERA5_VOR${lev[$nl]}_6hr_2000_DET_T42/tr_trs_pos

        utils/bin/box $file ${lats} ${latn} ${lonl} ${lonr} $option 0 0.0 
        mv ${file}.new ${file}.${option}_${lats}${latn}-${lonl}${lonr}
        utils/bin/tr2nc ${file}.${option}_${lats}${latn}-${lonl}${lonr} s utils/TR2NC/tr2nc.meta
    done
fi

if [ $pro == 4 ];then
    echo "=========== match low level with high level ==================="
    cd ~/TRACK-1.5.2/
    pwd
    if [ $filt == 1 ]; then 
        dir=${OUTDIR}filt${option}_${lats}${latn}-${lonl}${lonr}match${lev[$nl1]}_${lev[$nl2]}
        file1=${OUTDIR}${file}_${lev[$nl1]}_1980-2020.${option}_${lats}${latn}-${lonl}${lonr}
        file2=${OUTDIR}${file}_${lev[$nl2]}_1980-2020.${option}_${lats}${latn}-${lonl}${lonr}
    else
        dir=${OUTDIR}match${lev[$nl1]}_${lev[$nl2]}
        file1=${OUTDIR}${file}_${lev[$nl1]}_1980-2020
        file2=${OUTDIR}${file}_${lev[$nl2]}_1980-2020
    fi
    mkdir $dir 
    utils/bin/censemble2 $file1 $file2 0 100 10 1 0 0 0 0 s -1 1 5.0 0.1 > ${dir}/record
    
    mv ./match_ens* $dir
    mv ./str.dat $dir
    mv ./diff.dat $dir
    mv ./temp-dist.stat $dir

    #utils/bin/tr2nc ${dir}/match_ens1_yes.dat s utils/TR2NC/tr2nc.meta
    #utils/bin/tr2nc ${dir}/match_ens1_no.dat s utils/TR2NC/tr2nc.meta
    #utils/bin/tr2nc ${dir}/match_ens2_yes.dat s utils/TR2NC/tr2nc.meta
    #utils/bin/tr2nc ${dir}/match_ens2_no.dat s utils/TR2NC/tr2nc.meta

    rename .dat.nc .nc $dir/*
    rename match_ens1 ${lev[$nl1]}match${lev[$nl2]} $dir/*
    rename match_ens2 ${lev[$nl2]}match${lev[$nl1]} $dir/*
fi

