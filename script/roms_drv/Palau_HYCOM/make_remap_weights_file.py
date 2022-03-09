import matplotlib
matplotlib.use('Agg')
import pyroms
import pyroms_toolbox
import os
os.environ['PROJ_LIB'] = '/home/metctm1/array/soft/anaconda3/share/proj'

print("load the grid")
hycom_data = '/home/lzhenn/drv_field/hycom_subset/2021091500/hycom.njord.grid.exp930.nc' 
roms_hist = '/home/lzhenn/cooperate/data/case_study/coupled/2020060412/njord_his_d01.20200604.nc'
roms_grid = '/home/lzhenn/Njord/Projects/Njord/roms_swan_grid/roms_d01_lp0d1.nc'
srcgrd = pyroms_toolbox.Grid_HYCOM.get_nc_Grid_HYCOM(hycom_data,name='GLBy0.08')
dstgrd = pyroms.grid.get_ROMS_grid('njord',hist_file=roms_hist,grid_file=roms_grid)

print("make remap grid file for scrip")
pyroms_toolbox.Grid_HYCOM.make_remap_grid_file(srcgrd)
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='rho')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='u')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='v')

print("compute remap weights for t to rho")
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_GLBy0.08_t.nc'
grid2_file = 'remap_grid_njord_rho.nc'
interp_file1 = 'remap_weights_GLBy0.08_to_njord_bilinear_t_to_rho.nc'
interp_file2 = 'remap_weights_njord_to_GLBy0.08_bilinear_rho_to_t.nc'
map1_name = 'GLBy0.08 to njord Bilinear Mapping'
map2_name = 'njord to GLBy0.08 Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)


print("compute remap weights for t to u")
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_GLBy0.08_t.nc'
grid2_file = 'remap_grid_njord_u.nc'
interp_file1 = 'remap_weights_GLBy0.08_to_njord_bilinear_t_to_u.nc'
interp_file2 = 'remap_weights_njord_to_GLBy0.08_bilinear_u_to_t.nc'
map1_name = 'GLBy0.08 to njord Bilinear Mapping'
map2_name = 'njord to GLBy0.08 Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)


print("compute remap weights for t to v")
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_GLBy0.08_t.nc'
grid2_file = 'remap_grid_njord_v.nc'
interp_file1 = 'remap_weights_GLBy0.08_to_njord_bilinear_t_to_v.nc'
interp_file2 = 'remap_weights_njord_to_GLBy0.08_bilinear_v_to_t.nc'
map1_name = 'GLBy0.08 to njord Bilinear Mapping'
map2_name = 'njord to GLBy0.08 Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)

