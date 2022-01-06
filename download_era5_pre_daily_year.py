#!/usr/bin/env python

import cdsapi
from cdo import *
import numpy as np

cdo = Cdo()
c = cdsapi.Client()

day28 = ['01', '02', '03','04', '05', '06', '07', '08', '09',
         '10', '11', '12', '13', '14', '15', '16', '17', '18',
         '19', '20', '21', '22', '23', '24', '25', '26', '27',
         '28',]

day29 = ['01', '02', '03','04', '05', '06', '07', '08', '09',
         '10', '11', '12', '13', '14', '15', '16', '17', '18',
         '19', '20', '21', '22', '23', '24', '25', '26', '27',
         '28','29',]

day30 = ['01', '02', '03','04', '05', '06', '07', '08', '09',
         '10', '11', '12', '13', '14', '15', '16', '17', '18',
         '19', '20', '21', '22', '23', '24', '25', '26', '27',
         '28', '29', '30',]

day31 = ['01', '02', '03','04', '05', '06', '07', '08', '09',
         '10', '11', '12', '13', '14', '15', '16', '17', '18',
         '19', '20', '21', '22', '23', '24', '25', '26', '27',
         '28', '29', '30', '31',]

mon_day = {1:day31, 2:day28, 3:day31, 4:day30, 5:day31, 6:day30, 7:day31, 8:day31, 9:day30, 10:day31, 11:day30, 12:day31, }
mon_day_29 = {1:day31, 2:day29, 3:day31, 4:day30, 5:day31, 6:day30, 7:day31, 8:day31, 9:day30, 10:day31, 11:day30, 12:day31, }
#mon_day = {11:day30, 12:day31, }


year_mon = {1979:mon_day, 1980:mon_day_29, 1981:mon_day, 1982:mon_day, 1983:mon_day, 1984:mon_day_29, 1985:mon_day, 1986:mon_day, 1987:mon_day, 1988:mon_day_29,
            1989:mon_day, 1990:mon_day, 1991:mon_day, 1992:mon_day_29, 1993:mon_day, 1994:mon_day, 1995:mon_day, 1996:mon_day_29, 1997:mon_day, 1998:mon_day,
            1999:mon_day, 2000:mon_day_29, 2001:mon_day, 2002:mon_day, 2003:mon_day, 2004:mon_day_29, 2005:mon_day, 2006:mon_day, 2007:mon_day, 2008:mon_day_29,
            2009:mon_day, 2010:mon_day, 2011:mon_day, 2012:mon_day_29, 2013:mon_day, 2014:mon_day, 2015:mon_day, 2016:mon_day_29, 2017:mon_day, 2018:mon_day,}
#mon_day = {11:day30, 12:day31, }

#variables = {'zg':'geopotential', 'q':'specific_humidity'}
variables = {'omega':'vertical_velocity'}

for var in ('omega',):
    for year in range(2015,2016):
        for mon in range(7,13):
            c.retrieve(
                'reanalysis-era5-pressure-levels',
                {
                    'product_type': 'reanalysis',
                    'format': 'netcdf',
                    'variable': variables[var],
                    'pressure_level': [
                        '1', '5', '10',
                        '20', '30', '50',
                        '70', '100','150',
                        '200', '250', '300',
                        '400', '500', '600',
                        '700', '850', '925',
                        '1000',
                    ],
                    'year': str(year),
                    'month': [
                        np.str(mon),
                    ],
                    'day': year_mon[year][mon],
                    'time': [
                        '00:00', '01:00', '02:00',
                        '03:00', '04:00', '05:00',
                        '06:00', '07:00', '08:00',
                        '09:00', '10:00', '11:00',
                        '12:00', '13:00', '14:00',
                        '15:00', '16:00', '17:00',
                        '18:00', '19:00', '20:00',
                        '21:00', '22:00', '23:00',
                    ],
                    'grid': [1.0, 1.0
                    ],
        
                },
                var+'_'+str(year)+'_'+np.str(mon)+'.nc')
        
            cdo.daymean(input = '/mnt/e/downloads/'+var+'_'+str(year)+'_'+np.str(mon)+'.nc', output = '/mnt/e/ERA5/'+var+'/'+var+'.'+str(year)+'-'+np.str(mon).zfill(2)+'.daily.nc')
    for year in range(2016,2019):
        for mon in range(1,13):
            c.retrieve(
                'reanalysis-era5-pressure-levels',
                {
                    'product_type': 'reanalysis',
                    'format': 'netcdf',
                    'variable': variables[var],
                    'pressure_level': [
                        '1', '5', '10',
                        '20', '30', '50',
                        '70', '100','150',
                        '200', '250', '300',
                        '400', '500', '600',
                        '700', '850', '925',
                        '1000',
                    ],
                    'year': str(year),
                    'month': [
                        np.str(mon),
                    ],
                    'day': year_mon[year][mon],
                    'time': [
                        '00:00', '01:00', '02:00',
                        '03:00', '04:00', '05:00',
                        '06:00', '07:00', '08:00',
                        '09:00', '10:00', '11:00',
                        '12:00', '13:00', '14:00',
                        '15:00', '16:00', '17:00',
                        '18:00', '19:00', '20:00',
                        '21:00', '22:00', '23:00',
                    ],
                    'grid': [1.0, 1.0
                    ],
        
                },
                var+'_'+str(year)+'_'+np.str(mon)+'.nc')
        
            cdo.daymean(input = '/mnt/e/downloads/'+var+'_'+str(year)+'_'+np.str(mon)+'.nc', output = '/mnt/e/ERA5/'+var+'/'+var+'.'+str(year)+'-'+np.str(mon).zfill(2)+'.daily.nc')

