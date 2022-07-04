#!/bin/bash
#SBATCH --partition=long-serial
#SBATCH --output=/home/users/qd201969/uor_track/TMP/slurm-%j.out 
#SBATCH --error=/home/users/qd201969/uor_track/TMP/slurm-%j.err
#SBATCH --job-name=renql
#SBATCH --time=3-00:00:00

# combine
#sh ~/uor_track/control_era5_1hr_track.sh fft 1 0 0 25 45 60 110
#sh ~/uor_track/control_era5_1hr_track.sh ff 1 0 0 25 45 60 110

./2203-calc_maximum10mwind_mpool.py

