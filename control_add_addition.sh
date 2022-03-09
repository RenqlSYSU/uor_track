#!/bin/bash
#SBATCH --partition=cluster
#SBATCH --output=/home/users/qd201969/ERA5-1HR-lev/TMP/slurm-%j.out 
#SBATCH --error=/home/users/qd201969/ERA5-1HR-lev/TMP/slurm-%j.err
#SBATCH --job-name=ew
#SBATCH --time=10-00:00:00

lev=(250 500 850)
OUTDIR=/home/users/qd201969/ERA5-1HR-lev/

cd /home/users/qd201969/TRACK-1.5.2
pwd
echo "add additional fields"

for ny in {1980..2020};do
    for nl in {0..2};do
        addz=${OUTDIR}ERA5_VOR${lev[$nl]}_1hr_${ny}_DET/ff_trs_pos.addZ
        addwind=${OUTDIR}ERA5_VOR${lev[$nl]}_1hr_${ny}_DET/ff_trs_pos.addZ_addwind10m
        addprep=${OUTDIR}ERA5_VOR${lev[$nl]}_1hr_${ny}_DET/ff_trs_pos.addZ_addwind10m_precip

        if [ ! -f $addz ]; then
            sed -e "s/YYYY/${ny}/;s/LEV/${lev[$nl]}/" addZ.in > addZ
            bin/track.linux < addZ
            mv outdat/ff_trs.linux_addfld $addz 
        fi
        
        if [ ! -f $addwind ]; then
        if [ ! -f $addz ];then
            echo "${addz} does not exist"  
        else
            sed -e "s/YYYY/${ny}/;s/LEV/${lev[$nl]}/" addwind.in > addwind
            bin/track.linux < addwind
            mv outdat/ff_trs.linux_addfld $addwind 
        fi
        fi

        #if [ ! -f $addprep ]; then
        #if [ ! -f $addwind ]; then
        #    echo "${addwind} does not exist"  
        #else
        sed -e "s/YYYY/${ny}/;s/LEV/${lev[$nl]}/" addprecip.in > addprecip
        bin/track.linux < addprecip
        mv outdat/ff_trs.linux_addfld $addprep 
        #fi
        #fi

        if [ -f $addprep ]; then
            num1=$(awk 'END {print NR}' $addz)
            num2=$(awk 'END {print NR}' $addprep)
            if [ $num1 -eq $num2 ]; then
                echo "${addprep} right"
            else
                echo "${addprep} wrong"
            fi
        else
            echo "${addprep} still does not exit"
        fi
    done
done

