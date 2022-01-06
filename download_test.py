import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels-daily-means',
    {
        'format': 'netcdf',
        'product_type': 'daily_averaged_reanalysis',
        'variable': 'surface_pressure',
        'year': '1979',
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
    },
    '/home/users/qd201969/data/ERA5_test.nc')

