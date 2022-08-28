#!/bin/bash
#SBATCH --partition=par-single
#SBATCH --ntasks=3
#SBATCH --output=/home/users/qd201969/uor_track/TMP/slurm-%j.out 
#SBATCH --error=/home/users/qd201969/uor_track/TMP/slurm-%j.err
#SBATCH --job-name=renql
#SBATCH --time=2-00:00:00

#./2203-calc_maximum10mwind_mpool.py
#./2203-calc_clim_precip_mpool.py

# combine
#sh ~/uor_track/control_era5_1hr_track.sh fft 1 0 0 25 45 60 110
#sh ~/uor_track/control_era5_1hr_track.sh ff -1 0 0 25 45 60 110
sh ~/uor_track/control_era5_1hr_track.sh tr 2 0 '' 

#lev=(850 500 250)
#OUTDIR=/home/users/qd201969/ERA5-1HR-lev/
#for ny in {1980..2020};do
#    for nl in {0..2};do
#        rm ${OUTDIR}ERA5_VOR${lev[$nl]}_1hr_${ny}_DET/ff_trs_pos.addZ_addwind10m_precip-1211
#        rm ${OUTDIR}ERA5_VOR${lev[$nl]}_1hr_${ny}_DET/fft_trs_pos.addZ_addwind10m_precip-1211
#    done
#done
