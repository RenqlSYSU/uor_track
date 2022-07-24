#!/usr/bin/env python
import cdsapi

varname = ['Total precipitation']
filname = ['tp']
outdir = '/home/lzhenn/cooperate/data/cwrf_s2s'

c = cdsapi.Client()

for nv in range(0,len(varname),1):
    for ny in range(2011,2021):
        c.retrieve(
            'reanalysis-era5-single-levels-monthly-means',
            {
                'format': 'netcdf',
                'product_type': 'monthly_averaged_reanalysis',
                'variable': varname[nv],
                'year': str(ny),
                'month': '07',
                'time': '00:00',
                'expver': '1',
                'grid': [0.25, 0.25],
                'area': [50, 75, 15, 137, ],
            },
            '%s/ERA5-%s_%d.nc'%(outdir,filname[nv],ny))

