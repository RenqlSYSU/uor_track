#!/usr/bin/env python
import cdsapi, os

varname = ['u_component_of_wind','v_component_of_wind','temperature','geopotential']
filname = ['u','v','t','z']
path = '/gws/nopw/j04/ncas_generic/users/renql/ERA5_subdaily/'

c = cdsapi.Client()

for nv in range(0,len(varname),1):
    if not os.path.exists(path+filname[nv]):
        os.makedirs(path+filname[nv])

    for year in range(1979,2020,1):
        c.retrieve(
            'reanalysis-era5-pressure-levels',
            {
                'product_type': 'reanalysis',
                'expver': '1',
                'variable': varname[nv],
                'pressure_level': [
                    '200', '225', '250',
                    '450', '500', '825',
                    '850',
                ],
                'grid': [1.0, 1.0],
                'year': str(year), 
                'month': [
                    '01', '02', '03',
                    '04', '05', '06',
                    '07', '08', '09',
                    '10', '11', '12',
                ],
                'day': [
                    '01', '02', '03',
                    '04', '05', '06',
                    '07', '08', '09',
                    '10', '11', '12',
                    '13', '14', '15',
                    '16', '17', '18',
                    '19', '20', '21',
                    '22', '23', '24',
                    '25', '26', '27',
                    '28', '29', '30',
                    '31',
                ],
                'time': [
                    '00:00', '06:00', '12:00',
                    '18:00',
                ],
                'area': [
                    90, 0, 0,
                    150,
                ],
                'format': 'netcdf',
            },
            '%s%s/ERA5_NH_%s_%d.nc'%(path,filname[nv],filname[nv],year))

'''
'00:00', '01:00', '02:00',
'03:00', '04:00', '05:00',
'06:00', '07:00', '08:00',
'09:00', '10:00', '11:00',
'12:00', '13:00', '14:00',
'15:00', '16:00', '17:00',
'18:00', '19:00', '20:00',
'21:00', '22:00', '23:00',
'''