import xarray as xr
import xesmf

def regrid_GLBy(fld, method='bilinear'):
    hycom_data = '/home/lzhenn/drv_field/hycom_subset/2021091500/hycom_glby_930_2021091500.nc'
    roms_grid = '/home/lzhenn/Njord/Projects/Njord/roms_swan_grid/roms_d01_lp0d1.nc'
    coords = xr.open_dataset(roms_grid)
    coords = coords.rename({'lon_rho': 'lon', 'lat_rho': 'lat'})
    gsource = xr.open_dataset(hycom_data)

    regrid = xesmf.Regridder(
        gsource,
        coords,
        method=method,
    )
    tdest = regrid(fld)
    return tdest
        #filename='regrid_t.nc',
        #periodic=False,
        #reuse_weights=False
