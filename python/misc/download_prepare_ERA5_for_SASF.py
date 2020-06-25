#!/usr/bin/env python                                                                                                   
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-                                              

import sys
from os import path
import cdsapi
from netCDF4 import Dataset
from math import copysign
import numpy as nmp

# Coordinates (point) we want to extract for STATION ASF:

plon = 36.75 ; plat = 81. ; # East of Svalbard

list_crd_expected = ['longitude', 'latitude', 'time']
# Their name in the downloaded file:
list_var_expected = ['u10', 'v10', 'd2m', 'fal', 't2m', 'istl1', 'msl', 'siconc', 'sst', 'skt', 'ssrd', 'strd', 'tp']
# Their name in the cdsapi request:
lvdl = [ '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_dewpoint_temperature', 'forecast_albedo', '2m_temperature', 'ice_temperature_layer_1', 'mean_sea_level_pressure','sea_ice_cover', 'sea_surface_temperature', 'skin_temperature', 'surface_solar_radiation_downwards', 'surface_thermal_radiation_downwards', 'total_precipitation' ]
yyyy = 2018

# In output file:
cv_lon = 'nav_lon'
cv_lat = 'nav_lat'
cv_tim = 'time_counter'

# List of flux variables to convert to right unit (divide by rdt):
rdt = 3600.
list_flx = ['ssrd',       'strd',     'tp']
list_fnu = ['W m**-2', 'W m**-2',  'mm s**-1' ]
fact_flx = [ 1./rdt,    1./rdt ,   1000./rdt ] ; # tp in 'm' ...

list_temp_to_degC = [ 'sst' , 'skt' ] ; # For some reasons SAS part of NEMO expects SSTs in deg. Celsius...

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def long_to_m180_p180(xx):
    ## Forces longitude to be in the -180:180 frame...
    ## xx: longitude
    xx   = xx % 360.
    rlon = copysign(1.,180.-xx)*min(xx, abs(xx-360.)) ; 
    return rlon



nbfld = len(list_var_expected)
if len(lvdl) != nbfld:
    print(' ERROR: download list "lvdl" not the same size as "list_var_expected"!!!') ; sys.exit(0)

# Coordinates of 10deg-wide box to download:
dw = 5. ; # degrees
irng_lon = [ int(round(long_to_m180_p180(plon-dw),0)) , int(round(long_to_m180_p180(plon+dw),0)) ]
irng_lat = [ int(round(max(plat-dw,-90.),0))          , int(round(min(plat+dw, 90.),0))          ]

plon = long_to_m180_p180(plon)

print(' * Longitude =', plon, ' =>', irng_lon[:])
print(' * Latitude  =', plat, ' =>', irng_lat[:])
print('')




# Downloading, month after month...

c = cdsapi.Client()

for jm in range(12):

    cm = '%2.2i'%(jm+1)

    cf_fi = 'ERA5_arctic_surface_'+str(irng_lat[1])+'-'+str(irng_lat[0])+'_'+str(irng_lon[0])+'-'+str(irng_lon[1])+'_'+str(yyyy)+cm+'.nc'
    cf_fo = 'ERA5_arctic_surface_'+str(plat)+'N_'+str(plon)+'E_1h_'+str(yyyy)+cm+'.nc'
    
    if not path.exists(cf_fi):
    
        print('\nDoing month '+cm+' !')
    
        c.retrieve(
            'reanalysis-era5-single-levels',
            {
                'product_type': 'reanalysis',
                'format': 'netcdf',
                'variable': [
                    lvdl[0], lvdl[1], lvdl[2], lvdl[3], lvdl[4], lvdl[5], lvdl[6], lvdl[7], lvdl[8], lvdl[9], lvdl[10], lvdl[11], lvdl[12],
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
                    irng_lat[1], irng_lon[0], irng_lat[0],
                    irng_lon[1],
                ],
            },
            cf_fi )
    
    else:
        print('\nAlready done month '+cm+' !')

    print('')




    # Gonna fix this crap!
    
    id_fi = Dataset(cf_fi)

    # 1/ populate variables and check it's what's expected:
    list_var = id_fi.variables.keys()
    print(' *** list_var =', list_var)
    if list_var[:3] != list_crd_expected:
        print(' ERROR this is not the list of coordinates we expected...') ; sys.exit(0)
    if list_var[3:] != list_var_expected:
        print(' ERROR this is not the list of variables we expected...') ; sys.exit(0)

    Ni = id_fi.dimensions['longitude'].size
    Nj = id_fi.dimensions['latitude'].size
    #if not id_fi.dimensions['time'].isunlimited(): print 'PROBLEM: the time dimension is not UNLIMITED! Bad!'; sys.exit(0)
    Nt = id_fi.dimensions['time'].size ; # Not unlimited in downloaded files...
    print ' *** Input file: Ni, Nj, Nt = ', Ni, Nj, Nt, '\n'
        
    vlon  = id_fi.variables['longitude'][:] ; cunt_lon = id_fi.variables['longitude'].units ; clnm_lon = id_fi.variables['longitude'].long_name
    vlat  = id_fi.variables['latitude'][:]  ; cunt_lat = id_fi.variables['latitude'].units  ; clnm_lat = id_fi.variables['latitude'].long_name
    vtime = id_fi.variables['time'][:]      ; cunt_tim = id_fi.variables['time'].units      ; clnm_tim = id_fi.variables['time'].long_name

    ip = nmp.argmin(nmp.abs(vlon-plon))
    jp = nmp.argmin(nmp.abs(vlat-plat))

    print(' *** ip, jp =', ip, jp)


    # Creating output file for ocean:
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ni = 3 ; nj = 3
    id_fo = Dataset(cf_fo, 'w', format='NETCDF4')

    # Dimensions:
    id_fo.createDimension('x'   , ni  )
    id_fo.createDimension('y'   , nj  )
    id_fo.createDimension(cv_tim, None)

    # Variables
    ido_lon = id_fo.createVariable(cv_lon, 'f4', ('y','x',), zlib=True) ; ido_lon.units = cunt_lon ; ido_lon.long_name = clnm_lon
    ido_lat = id_fo.createVariable(cv_lat, 'f4', ('y','x',), zlib=True) ; ido_lat.units = cunt_lat ; ido_lat.long_name = clnm_lat
    ido_tim = id_fo.createVariable(cv_tim, 'f4', (cv_tim,) , zlib=True) ; ido_tim.units = cunt_tim ; ido_tim.long_name = clnm_tim
    
    # Creating fields in output file:
    ido_var = []
    iv = 0
    for cvar in list_var_expected:        
        ido_var.append(id_fo.createVariable(cvar, 'f4', (cv_tim,'y','x',), zlib=True))
        if cvar in list_flx:
            idx = list_flx.index(cvar)
            ido_var[iv].units     = list_fnu[idx]
        elif cvar in list_temp_to_degC:
            ido_var[iv].units     = 'degC'        
        else:
            ido_var[iv].units = id_fi.variables[cvar].units
        ido_var[iv].long_name = id_fi.variables[cvar].long_name
        iv = iv + 1
    
    # Filling coordinates:
    ido_lon[:,:] = vlon[ip]
    ido_lat[:,:] = vlat[jp]

    # Filling fields
    for jt in range(Nt):
        ido_tim[jt] = vtime[jt]
        iv = 0
        for cvar in list_var_expected:
            ido_var[iv][jt,:,:] = id_fi.variables[cvar][jt,jp,ip]
            #
            # Flux conversion ???
            if cvar in list_flx:
                idx = list_flx.index(cvar)
                ido_var[iv][jt,:,:] = ido_var[iv][jt,:,:] * fact_flx[idx]
            if cvar in list_temp_to_degC:
                ido_var[iv][jt,:,:] = ido_var[iv][jt,:,:] - 273.15

            #
            iv = iv + 1


    id_fo.About = "Input file for 'STATION_ASF' NEMO test-case, generated with 'download_prepare_ERA5_for_SASF.py' of AeroBulk (https://github.com/brodeau/aerobulk)."
            
    id_fi.close()
    id_fo.close()
