#!/usr/bin/env python3                                                                                                   
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-                                              
#
# LOLO: add an empty array for SSH in output netcdf file !!!
##      + same for ssu and ssv !
##      + add constant salinity sss of say 34 PSU
##

import sys
from os import path
import cdsapi
from netCDF4 import Dataset,num2date
from math import copysign
import numpy as nmp

yyyy = 2018

l_interp_daily = True

# Coordinates (point) we want to extract for STATION ASF:

plon = -36.  ; plat = 84. ; # North Greenland...
#plon = 36.75 ; plat = 81. ; # East of Svalbard
#plon = -65.1  ; plat = 73.2 ; # Center of Baffin Bay

list_crd_expected = ['longitude', 'latitude', 'time']
# Their name in the downloaded file:
### dumping 'fal','forecast_albedo' since it's proportional to ice fraction...
list_var_expected = ['u10', 'v10', 'd2m', 't2m', 'istl1', \
                     'msl', 'skt','ssrd', 'strd', 'tp', 'sf' ]

# Their name in the cdsapi request:
lvdl_h = [ '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_dewpoint_temperature', '2m_temperature', \
           'ice_temperature_layer_1', 'mean_sea_level_pressure', 'skin_temperature', 'surface_solar_radiation_downwards', \
           'surface_thermal_radiation_downwards', 'total_precipitation', 'snowfall' ]

# Daily fields:
list_var_daily_expected = [    'siconc'    ,        'sst'              ];#, ] 
lvdl_d                  = [ 'sea_ice_cover', 'sea_surface_temperature' ];#, ] 



# In output file:
cv_lon = 'nav_lon'
cv_lat = 'nav_lat'
cv_tim = 'time_counter'

# List of flux variables to convert to right unit (divide by rdt):
rdt = 3600.
list_flx = ['ssrd',       'strd',     'tp'   ,     'sf'     ]
list_fnu = ['W m**-2', 'W m**-2',  'mm s**-1',  'mm s**-1'  ]
fact_flx = [ 1./rdt,    1./rdt  ,  1000./rdt ,  1000./rdt   ] ; # tp and sf in 'm' ...

list_temp_to_degC = [ 'sst' , 'skt' ] ; # For some reasons SAS part of NEMO expects SSTs in deg. Celsius...

# Extra fields (for the SAS part) to add and their value (constant along time records...):
list_extra = [ 'sss'                 , 'ssh'               ,  'ssu'                 , 'ssv'                       ,    'ialb'        ]
rval_extra = [  34.                  ,   0.                ,     0.                 ,   0.                        ,      0.55        ]
cunt_extra = [  ''                   ,  'm'                ,    'm s**-1'           , 'm s**-1'                   ,       ''         ]
clnm_extra = [ 'Sea surface salinity', 'Sea surface height', 'Zonal surface current', 'Meridional surface current', 'Sea-ice albedo' ]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def long_to_m180_p180(xx):
    ## Forces longitude to be in the -180:180 frame...
    ## xx: longitude
    xx   = xx % 360.
    rlon = copysign(1.,180.-xx)*min(xx, abs(xx-360.)) ; 
    return rlon


# Coordinates of 10deg-wide box to download:
dw = 5. ; # degrees
irng_lon = [ int(round(long_to_m180_p180(plon-dw),0)) , int(round(long_to_m180_p180(plon+dw),0)) ]

irng_lat = [ int(round(max(plat-dw,-90.),0))          , int(round(min(plat+dw, 90.),0))          ]

plon = long_to_m180_p180(plon)

print(' * Longitude =', plon, ' =>', irng_lon[:])
print(' * Latitude  =', plat, ' =>', irng_lat[:])
print('')
cinfo_box5x5 = str(irng_lat[1])+'-'+str(irng_lat[0])+'_'+str(irng_lon[0])+'-'+str(irng_lon[1])
cinfo_coord  = str(plat)+'N_'+str(plon)+'E'






# Monthly stuff:
nbfld_daily = len(list_var_daily_expected)
if len(lvdl_d) != nbfld_daily:
    print(' ERROR: download list "lvdl_d" not the same size as "list_var_daily_expected"!!!', len(lvdl_d), nbfld_daily) ; sys.exit(0)

# Hourly stuff:
nbfld = len(list_var_expected)
if len(lvdl_h) != nbfld:
    print(' ERROR: download list "lvdl_h" not the same size as "list_var_expected"!!!', len(lvdl_h), nbfld)
    print(lvdl_h,'\n')
    print(list_var_expected,'\n')
    sys.exit(0)

nbfld_tot = nbfld
if l_interp_daily: nbfld_tot = nbfld + nbfld_daily



    
c = cdsapi.Client()


#################
# Daily fields #
#################

Nt_tot_daily = 0

for jm in range(12):

    cm = '%2.2i'%(jm+1)

    cf_fi_d = 'ERA5_arctic_surface__BOX-5x5deg_'+cinfo_box5x5+'__'+str(yyyy)+cm+'_daily.nc'
    #cf_fo_d = 'ERA5_arctic_surface_'+cinfo_coord+'_1d_y'+str(yyyy)+'m'+cm+'.nc'

    print(' *** cf_fi_d = '+cf_fi_d)
    #print(' *** cf_fo_d = '+cf_fo_d)
    
    if not path.exists(cf_fi_d):
    
        print('\nDoing month '+cm+' !')
    
        c.retrieve(
            'reanalysis-era5-single-levels',
            {
                'product_type': 'reanalysis',
                'format': 'netcdf',
                'variable': lvdl_d,
                'year': str(yyyy),
                'month': [ cm, ],
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
                'time': [ '12:00' ],
                'area': [
                    irng_lat[1], irng_lon[0], irng_lat[0],
                    irng_lon[1], # 
                ],
            },
            cf_fi_d )
        
        if cm == '01':
            c.retrieve(
                'reanalysis-era5-single-levels',
                {
                    'product_type': 'reanalysis',
                    'format': 'netcdf',
                    'variable': lvdl_d,
                    'year': str(yyyy-1),
                    'month': [ '12', ],
                    'day': [ '31', ],
                    'time': [ '12:00' ],
                    'area': [
                        irng_lat[1], irng_lon[0], irng_lat[0],
                        irng_lon[1],
                    ],
                },
                cf_fi_d+'.before' )            
            
        if cm == '12':
            c.retrieve(
                'reanalysis-era5-single-levels',
                {
                    'product_type': 'reanalysis',
                    'format': 'netcdf',
                    'variable': lvdl_d,
                    'year': str(yyyy+1),
                    'month': [ '01', ],
                    'day': [ '31', ],
                    'time': [ '12:00' ],
                    'area': [
                        irng_lat[1], irng_lon[0], irng_lat[0],
                        irng_lon[1],
                    ],
                },
                cf_fi_d+'.after' )
            
        
    else:
        print('\nAlready done month '+cm+' !')

    
    print('')

    # Gonna fix this crap! #lulu

    list_fi_check = [ cf_fi_d ]
    if cm == '01': list_fi_check = [ cf_fi_d+'.before', cf_fi_d ]
    if cm == '12': list_fi_check = [ cf_fi_d, cf_fi_d+'.after' ]
    
    for ff in list_fi_check:
        print('\n Checking file '+ff+' !')
        id_fi = Dataset(ff)
        # 1/ populate variables and check it's what's expected:
        list_var = list(id_fi.variables.keys())
        print(' *** list_var          =', list_var)    
        if list_var[:3] != list_crd_expected:
            print(' ERROR this is not the list of coordinates we expected...') ; sys.exit(0)
        if list_var[3:] != list_var_daily_expected:
            print(' ERROR this is not the list of variables we expected...') ; sys.exit(0)
    
        Ni = id_fi.dimensions['longitude'].size
        Nj = id_fi.dimensions['latitude'].size
        #if not id_fi.dimensions['time'].isunlimited(): print 'PROBLEM: the time dimension is not UNLIMITED! Bad!'; sys.exit(0)
        Nt = id_fi.dimensions['time'].size ; # Not unlimited in downloaded files...
        print(' *** Input file: Ni, Nj, Nt = ', Ni, Nj, Nt, '\n')
        Nt_tot_daily = Nt_tot_daily + Nt
        #vtime_d = id_fi.variables['time'][:]
        
        id_fi.close()



# The 12 month for the daily field have been downloaded
print('\n Number of days for year '+str(yyyy)+': '+str(Nt_tot_daily-2))

xdata_d = nmp.zeros((nbfld_daily,Nt_tot_daily))  ; # before and after
vtime_d = nmp.zeros(             Nt_tot_daily )



# Reading the 12 files and filling the vtime_d and xdata_d (whole time-series for current year):
jt = 0
for jm in range(12):
    cm = '%2.2i'%(jm+1)
    cf_fi_d = 'ERA5_arctic_surface__BOX-5x5deg_'+cinfo_box5x5+'__'+str(yyyy)+cm+'_daily.nc'

    if cm == '01':
        id_fi = Dataset(cf_fi_d+'.before')
        vt = id_fi.variables['time'][:]
        if len(vt) != 1: print('ERROR #1!'); sys.exit(1)
        if jt==0: cunit_t = id_fi.variables['time'].units
        #
        vlon  = id_fi.variables['longitude'][:]
        vlat  = id_fi.variables['latitude'][:]
        ip = nmp.argmin(nmp.abs(vlon-plon))
        jp = nmp.argmin(nmp.abs(vlat-plat))
        print(' *** ip, jp =', ip, jp)
        #
        vtime_d[jt] = vt[0]
        jv = 0
        for cv in list_var_daily_expected:
            xdata_d[jv,jt] = id_fi.variables[cv][0,jp,ip]
            jv=jv+1
        id_fi.close()
        jt = jt + 1

    # Always:
    id_fi = Dataset(cf_fi_d)
    vt = id_fi.variables['time'][:]
    Nt = len(vt)
    vtime_d[jt:Nt+jt] = vt[:]
    jv = 0
    for cv in list_var_daily_expected:
        xdata_d[jv,jt:Nt+jt] = id_fi.variables[cv][:,jp,ip]
        jv=jv+1
    id_fi.close()
    jt = jt + Nt

    if cm == '12':
        id_fi = Dataset(cf_fi_d+'.after')
        vt = id_fi.variables['time'][:]
        if len(vt) != 1: print('ERROR #1!'); sys.exit(2)
        vtime_d[jt] = vt[0]
        jv = 0
        for cv in list_var_daily_expected:
            xdata_d[jv,jt] = id_fi.variables[cv][0,jp,ip]
            jv=jv+1
        id_fi.close()
        jt = jt + 1


print('\n jt, Nt_tot_daily =', jt, Nt_tot_daily)

# Debug to control...
print('\n\n')
for jt in range(Nt_tot_daily):
    print(jt, vtime_d[jt], num2date(vtime_d[jt], units=cunit_t))
print('\n\n')


## ===> so vtime_d[:] and xdata_d[:,:] is what needs to be interpolated later !


#sys.exit(0)







#################
# Hourly fields #
#################
    
for jm in range(12):

    cm = '%2.2i'%(jm+1)
    
    cf_fi = 'ERA5_arctic_surface__BOX-5x5deg_'+cinfo_box5x5+'__'+str(yyyy)+cm+'.nc'
    cf_fo = 'ERA5_arctic_surface_'+cinfo_coord+'_1h_y'+str(yyyy)+'m'+cm+'.nc'
    
    if not path.exists(cf_fi):
        
        print('\nDoing month '+cm+' !')
        
        c.retrieve(
            'reanalysis-era5-single-levels',
            {
                'product_type': 'reanalysis',
                'format': 'netcdf',
                'variable': lvdl_h,
                'year': str(yyyy),
                'month': [ cm, ],
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
    list_var = list(id_fi.variables.keys())
    print(' *** list_var =', list_var)
    if list_var[:3] != list_crd_expected:
        print(' ERROR this is not the list of coordinates we expected...') ; sys.exit(0)
    if list_var[3:] != list_var_expected:
        print(' ERROR this is not the list of variables we expected...') ; sys.exit(0)

    Ni = id_fi.dimensions['longitude'].size
    Nj = id_fi.dimensions['latitude'].size
    #if not id_fi.dimensions['time'].isunlimited(): print 'PROBLEM: the time dimension is not UNLIMITED! Bad!'; sys.exit(0)
    Nt = id_fi.dimensions['time'].size ; # Not unlimited in downloaded files...
    print(' *** Input file: Ni, Nj, Nt = ', Ni, Nj, Nt, '\n')
        
    vlon  = id_fi.variables['longitude'][:] ; cunt_lon = id_fi.variables['longitude'].units ; clnm_lon = id_fi.variables['longitude'].long_name
    vlat  = id_fi.variables['latitude'][:]  ; cunt_lat = id_fi.variables['latitude'].units  ; clnm_lat = id_fi.variables['latitude'].long_name
    vtime = id_fi.variables['time'][:]      ; cunt_tim = id_fi.variables['time'].units      ; clnm_tim = id_fi.variables['time'].long_name

    ip = nmp.argmin(nmp.abs(vlon-plon))
    jp = nmp.argmin(nmp.abs(vlat-plat))

    print(' *** ip, jp =', ip, jp)


    # Creating output file for ocean:
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ni = 3 ; nj = 3
    jr1 = 0 ; jr2 = 0
    
    id_fo = Dataset(cf_fo, 'w', format='NETCDF4')

    # Dimensions:
    id_fo.createDimension('x'   , ni  )
    id_fo.createDimension('y'   , nj  )
    id_fo.createDimension(cv_tim, None)

    # Variables
    ido_lon = id_fo.createVariable(cv_lon, 'f4', ('y','x',), zlib=True) ; ido_lon.units = cunt_lon ; ido_lon.long_name = clnm_lon # 
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

    # Dayly fields that will be interpolated in here:
    if l_interp_daily: #lulu
        for cvar in list_var_daily_expected:
            ido_var.append(id_fo.createVariable(cvar, 'f4', (cv_tim,'y','x',), zlib=True))
            if cvar in list_temp_to_degC:
                ido_var[iv].units = 'degC'
            else:
                ido_var[iv].units = 'boo'
            ido_var[iv].long_name = 'boo'
            iv = iv + 1
        
    # Creating extra fields:
    for cvar in list_extra:
        ido_var.append(id_fo.createVariable(cvar, 'f4', (cv_tim,'y','x',), zlib=True))
        ido_var[iv].units = cunt_extra[iv-nbfld_tot]
        ido_var[iv].long_name = clnm_extra[iv-nbfld_tot]
        iv = iv + 1
        
    # Filling coordinates:
    ido_lon[:,:] = vlon[ip]
    ido_lat[:,:] = vlat[jp]

    # Filling fields
    for jt in range(Nt):

        rt = vtime[jt] ; # current time!        
        ido_tim[jt] = rt
        
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

        # Dayly fields that will be interpolated in here:
        if l_interp_daily:
            #
            #print('LOLO: current time is rt =', rt)
            for jtd in range(jr1, Nt_tot_daily-1):
                if vtime_d[jtd] <= rt and vtime_d[jtd+1] > rt:
                    jr1 = jtd ; jr2 = jtd + 1 
                    break
            #print('LOLO: what we found is:', vtime_d[jr1], vtime_d[jr2])
            #
            ivd = 0
            for cvar in list_var_daily_expected:
                # Linear interpolation !!!
                rslope = (xdata_d[ivd,jr2] - xdata_d[ivd,jr1])/(vtime_d[jr2] - vtime_d[jr1])
                ido_var[iv][jt,:,:] = xdata_d[ivd,jr1] +  rslope*(rt - vtime_d[jr1])
                if cvar in list_temp_to_degC: ido_var[iv][jt,:,:] = ido_var[iv][jt,:,:] - 273.15
                iv = iv + 1 ; ivd = ivd + 1
                
        for cvar in list_extra:
            ido_var[iv][jt,:,:] = rval_extra[iv-nbfld_tot]
            iv = iv + 1

            

    id_fo.About = "Input file for 'STATION_ASF' NEMO test-case, generated with 'download_prepare_ERA5_for_SASF.py' of AeroBulk (https://github.com/brodeau/aerobulk)."
            
    id_fi.close()
    id_fo.close()
