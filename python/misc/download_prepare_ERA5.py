#!/usr/bin/env python3
#
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##
# Download ERA5 on the horizontal domain of your choice at 1h or 3h frequency...
#   Surface atmospheric variables required for atmospheric forcing of an OGCM...
#
## L. Brodeau, 2020
##
#
import sys
from os import path
import cdsapi
from netCDF4 import Dataset
from math import copysign
import numpy as nmp


dir_out = '/mnt/meom/workdir/brodeau/FATM/ERA5' ; # where to save...

cfrq="1h"
#cfrq="3h" ! => achtung for the multiplication for fluxes, this cannot be trusted!!!

l_download_only = False

# Coordinate box that contains your regional box of interest:
plat_min = -50. ; plon_min=  140.
plat_max =  35. ; plon_max=  -69.


list_crd_expected = ['longitude', 'latitude', 'time']
# Their name in the downloaded file:
list_var_expected = ['u10', 'v10', 'd2m', 't2m', 'msl', 'ssrd', 'strd', 'tp' ]
# Their name in the cdsapi request:
lvdl = [ '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_dewpoint_temperature', '2m_temperature', 'mean_sea_level_pressure', 'surface_solar_radiation_downwards', 'surface_thermal_radiation_downwards', 'total_precipitation' ]

vd_months = [ '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', ]

vd_days   = [ '01', '02', '03', '04', '05', '06', '07', '08', '09', '10',
              '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
              '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', ]

vd_snaps  =[ '00:00', '01:00', '02:00',
             '03:00', '04:00', '05:00',
             '06:00', '07:00', '08:00',
             '09:00', '10:00', '11:00',
             '12:00', '13:00', '14:00',
             '15:00', '16:00', '17:00',
             '18:00', '19:00', '20:00',
             '21:00', '22:00', '23:00', ]

rdt = 3600.


if cfrq  == "1h":
    
    # Hourly setup:
    # ~~~~~~~~~~~~
    vd_snaps  =[ '00:00', '01:00', '02:00', '03:00', '04:00', '05:00',
                 '06:00', '07:00', '08:00', '09:00', '10:00', '11:00',
                 '12:00', '13:00', '14:00', '15:00', '16:00', '17:00',
                 '18:00', '19:00', '20:00', '21:00', '22:00', '23:00', ]
    rdt = 3600. ; # for conversion from accumulated fluxes to fluxes (ex: W/m^2*s -> W/m^2)
    
elif cfrq  == "3h":
    
    # 3-Hourly setup:
    # ~~~~~~~~~~~~
    
    vd_snaps  =[ '00:00', '03:00',
                 '06:00', '09:00',
                 '12:00', '15:00',
                 '18:00', '21:00', ]
    rdt = 3600. ; # for conversion from accumulated fluxes to fluxes (ex: W/m^2*s -> W/m^2)
    #             # WARNING:  yeah I know it's weird it shoud be "3*3600"... but no!
    
else:
    print('\n UNKNOWN frequency!!! => '+cfrq+'\n')
    sys.exit(0)
    


def usage():
    print('\nUsage: '+sys.argv[0]+' <year> (<variable to treat>)')
    print('         *** known variables:'); print(list_var_expected[:],'\n')
    return sys.exit(0)



narg = len(sys.argv)
if not narg in [2,3]:
    usage()
yyyy = int(sys.argv[1])
#
if narg == 3:
    cvec = sys.argv[2]
    if not cvec in list_var_expected: usage()
    jv = list_var_expected.index(cvec)
    cvsc = lvdl[jv]
    # Updating lists with only the selected variable!
    del list_var_expected, lvdl

    list_var_expected = [ cvec ]
    lvdl              = [ cvsc ]
    print('   ===> gonna treat: ', list_var_expected[0], lvdl[0], '\n')

    
# In output file:
cv_lon = 'lon'
cv_lat = 'lat'
cv_tim = 'time'

# List of flux variables to convert to right unit (divide by rdt):
list_flx = ['ssrd',       'strd',     'tp'    ]
list_fnu = ['W m**-2', 'W m**-2',  'mm s**-1' ]
fact_flx = [ 1./rdt,    1./rdt  ,  1000./rdt  ] ; # tp and sf in 'm' ...





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def long_to_m180_p180(xx):
    ## Forces longitude to be in the -180:180 frame...
    ## xx: longitude
    xx   = xx % 360.
    rlon = copysign(1.,180.-xx)*min(xx, abs(xx-360.)) ;
    return rlon


nbfld = len(list_var_expected)
print('\n *** Going to download '+str(nbfld)+' fields!\n')

if len(lvdl) != nbfld:
    print(' ERROR: download list "lvdl" not the same size as "list_var_expected"!!!', len(lvdl), nbfld) ; sys.exit(0)

if (plon_min,plon_max) != (-180.,180.):
    plon_min = long_to_m180_p180(plon_min) ; plon_max = long_to_m180_p180(plon_max)

clab_coor = '_'+str(int(plat_min))+'N_'+str(int(plon_min))+'E_'+str(int(plat_max))+'N_'+str(int(plon_max))+'E_'

print('')
print(' * Longitude =', plon_min, plon_max)
print(' * Latitude  =', plat_min, plat_max)
print(' *  label    => ', clab_coor,'\n')


# Downloading, month after month...

c = cdsapi.Client()

for jv in range(nbfld):

    cvsc = lvdl[jv] ;              # variable name in script call
    cvec = list_var_expected[jv] ; # variable name in netcdf file

    cf_fi = dir_out+'/'+cvec+'_ERA5_surface'+clab_coor+str(yyyy)+'.nc'
    cf_fo = dir_out+'/'+cvec+'_ERA5_tropico-box_surface_'+cfrq+'_y'+str(yyyy)+'.nc'

    print(' cf_fi = ', cf_fi)
    if not path.exists(cf_fi):

        print('\nWill download field '+cvec+' !')

        c.retrieve(
            'reanalysis-era5-single-levels',
            {
                'product_type': 'reanalysis',
                'format': 'netcdf',
                'variable': [ cvsc, ],
                'year': str(yyyy),
                'month': vd_months,
                'day':   vd_days,
                'time': vd_snaps,
                'area': [ plat_max, plon_min, plat_min, plon_max, ],
            },
            cf_fi )

    else:
        print('\nAlready dowloaded field '+cvec+' !\n\n')








    if not l_download_only and not path.exists(cf_fo):

        # Gonna fix this piecec of crap!

        id_fi = Dataset(cf_fi)

        Ni = id_fi.dimensions['longitude'].size
        Nj = id_fi.dimensions['latitude'].size
        Nt = id_fi.dimensions['time'].size ; # Not unlimited in downloaded files...
        print(' *** Input file for '+cvec+': Ni, Nj, Nt = ', Ni, Nj, Nt, '\n')


        vlon  = id_fi.variables['longitude'][:]
        cunt_lon = id_fi.variables['longitude'].units ; clnm_lon = id_fi.variables['longitude'].long_name

        vlat  = id_fi.variables['latitude'][:]
        cunt_lat = id_fi.variables['latitude'].units  ; clnm_lat = id_fi.variables['latitude'].long_name

        vtime = id_fi.variables['time'][:]
        cunt_tim = id_fi.variables['time'].units      ; clnm_tim = id_fi.variables['time'].long_name


        # Creating output file
        # ~~~~~~~~~~~~~~~~~~~~~

        ni = Ni
        nj = Nj

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        id_fo = Dataset(cf_fo, 'w', format='NETCDF4')

        # Dimensions:
        id_fo.createDimension('x'   , ni  )
        id_fo.createDimension('y'   , nj  )
        id_fo.createDimension(cv_tim, None)

        # Variables
        ido_lon = id_fo.createVariable(cv_lon, 'f4', ('x',), zlib=True, complevel=4) ; ido_lon.units = cunt_lon ; ido_lon.long_name = clnm_lon
        ido_lat = id_fo.createVariable(cv_lat, 'f4', ('y',), zlib=True, complevel=4) ; ido_lat.units = cunt_lat ; ido_lat.long_name = clnm_lat
        ido_tim = id_fo.createVariable(cv_tim, 'f4', (cv_tim,) , zlib=True, complevel=4) ; ido_tim.units = cunt_tim ; ido_tim.long_name = clnm_tim

        # Creating fields in output file:
        ido_var = id_fo.createVariable(cvec, 'f4', (cv_tim,'y','x',), zlib=True, complevel=4)
        if cvec in list_flx:
            idx = list_flx.index(cvec)
            ido_var.units     = list_fnu[idx]
        else:
            ido_var.units = id_fi.variables[cvec].units
            ido_var.long_name = id_fi.variables[cvec].long_name

        # Filling coordinates:
        ido_lon[:] = vlon[:]
        ido_lat[:] = nmp.flipud(vlat[:])

        # Filling fields
        for jt in range(Nt):

            print('   * "'+cvec+'" => jt = ', jt, '  /', Nt)
            ido_tim[jt] = vtime[jt]

            if jt==0: Xvar = nmp.zeros((nj,ni))

            Xvar[:,:] = id_fi.variables[cvec][jt,:,:]

            Xvar[:,:] = nmp.flipud(Xvar[:,:])

            ido_var[jt,:,:] = Xvar[:,:]

            # Flux conversion ???
            if cvec in list_flx:
                idx = list_flx.index(cvec)
                ido_var[jt,:,:] = ido_var[jt,:,:] * fact_flx[idx]
                #
        id_fo.About = "Generated with 'download_prepare_ERA5.py' of AeroBulk (https://github.com/brodeau/aerobulk)."

        id_fo.close()
        id_fi.close()
