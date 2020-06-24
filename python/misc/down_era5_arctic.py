#!/usr/bin/env python                                                                                                   
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-                                              

from os import path
import cdsapi

yyyy = 2018

c = cdsapi.Client()

for jm in range(12):

    cm = '%2.2i'%(jm+1)

    cf_out = 'ERA5_arctic_surface_'+str(yyyy)+cm+'.nc'

    if not path.exists(cf_out):
    
        print('\nDoing month '+cm+' !')
    
        c.retrieve(
            'reanalysis-era5-single-levels',
            {
                'product_type': 'reanalysis',
                'format': 'netcdf',
                'variable': [
                    '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_dewpoint_temperature', 'forecast_albedo', '2m_temperature', 'ice_temperature_layer_1', 'mean_sea_level_pressure','sea_ice_cover', 'sea_surface_temperature', 'skin_temperature', 'surface_solar_radiation_downwards', 'surface_thermal_radiation_downwards', 'total_precipitation',
                ],
                'year': str(yyyy),
                'month': [
                    cm,
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
                    '00:00', '01:00', '02:00',
                    '03:00', '04:00', '05:00',
                    '06:00', '07:00', '08:00',
                    '09:00', '10:00', '11:00',
                    '12:00', '13:00', '14:00',
                    '15:00', '16:00', '17:00',
                    '18:00', '19:00', '20:00',
                    '21:00', '22:00', '23:00',
                ],
                'area': [
                    85, -180, 75,
                    180,
                ],
            },
            cf_out )
    
    else:
        print('\nAlready done month '+cm+' !')

    print('')
