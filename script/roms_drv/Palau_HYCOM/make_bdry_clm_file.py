import matplotlib
matplotlib.use('Agg')

import subprocess
import os
import sys
from multiprocessing import Pool
from functools import partial
import numpy as np
#import pdb

#increase the maximum number of open files allowed
#import resource
#resource.setrlimit(resource.RLIMIT_NOFILE, (3000,-1))

import pyroms
import pyroms_toolbox

from remap_bdry import remap_bdry
from remap_bdry_uv import remap_bdry_uv
from remap_clm import remap_clm
from remap_clm_uv import remap_clm_uv

def interp_clm(file, src_grd, dst_grd):
    zeta = remap_clm(file, 'ssh', src_grd, dst_grd, dst_dir=dst_dir)
    dst_grd = pyroms.grid.get_ROMS_grid('njord', zeta=zeta)
    remap_clm(file, 'temp', src_grd, dst_grd, dst_dir=dst_dir)
    remap_clm(file, 'salt', src_grd, dst_grd, dst_dir=dst_dir)
    remap_clm_uv(file, src_grd, dst_grd, dst_dir=dst_dir)

    # merge file
    clim_file = dst_dir + file.rsplit('/')[-1][:-3] + '_clim_' + dst_grd.name + '.nc'

    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_ssh_clim_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-O', out_file, clim_file) 
    subprocess.call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_temp_clim_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, clim_file) 
    subprocess.call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_salt_clim_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, clim_file) 
    subprocess.call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_u_clim_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, clim_file) 
    subprocess.call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_v_clim_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, clim_file) 
    subprocess.call(command)

def interp_bdry(file, src_grd, dst_grd):
    zeta = remap_bdry(file, 'ssh', src_grd, dst_grd, dst_dir=dst_dir)
    dst_grd = pyroms.grid.get_ROMS_grid('njord', zeta=zeta)
    remap_bdry(file, 'temp', src_grd, dst_grd, dst_dir=dst_dir)
    remap_bdry(file, 'salt', src_grd, dst_grd, dst_dir=dst_dir)
    remap_bdry_uv(file, src_grd, dst_grd, dst_dir=dst_dir)

    # merge file
    bdry_file = dst_dir + file.rsplit('/')[-1][:-3] + '_bdry_' + dst_grd.name + '.nc'

    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_ssh_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-O', out_file, bdry_file) 
    subprocess.call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_temp_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, bdry_file) 
    subprocess.call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_salt_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, bdry_file) 
    subprocess.call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_u_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, bdry_file) 
    subprocess.call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_v_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, bdry_file) 
    subprocess.call(command)
    os.remove(out_file)

data_dir = '/home/lzhenn/drv_field/hycom_subset/2021091500/'
dst_dir='/home/lzhenn/drv_field/icbc/2021091500'

#year = int(sys.argv[1])
#lst_year = sys.argv[1:]
year = '2021'
lst_year = [year]
lst_file = []

for year in lst_year:
    year = np.str(year)
    command = 'ls ' + data_dir + 'hycom_glby_930_' + year + '*'
    lst = subprocess.check_output(command, shell=True)
    lst = lst.split()
    lst_file = lst_file + lst

print('Build OBC file from the following file list:')
print(lst_file)
print(' ')

hycom_data = '/home/lzhenn/drv_field/hycom_subset/2021091500/hycom.njord.grid.exp930.nc' 
roms_hist = '/home/lzhenn/cooperate/data/case_study/coupled/2020060412/njord_his_d01.20200604.nc'
roms_grid = '/home/lzhenn/Njord/Projects/Njord/roms_swan_grid/roms_d01_lp0d1.nc'
src_grd = pyroms_toolbox.Grid_HYCOM.get_nc_Grid_HYCOM(hycom_data,name='GLBy0.08')
dst_grd = pyroms.grid.get_ROMS_grid('njord',hist_file=roms_hist,grid_file=roms_grid)

processes = 4
p = Pool(processes)
# Trick to pass more than one arg
partial_interp_bdry = partial(interp_bdry, src_grd=src_grd, dst_grd=dst_grd)
results = p.map(partial_interp_bdry, lst_file)

partial_interp_clm = partial(interp_clm, src_grd=src_grd, dst_grd=dst_grd)
results = p.map(partial_interp_clm, lst_file)

