#!/usr/bin/env python
#
# L. Brodeau, June 2019

import sys
import numpy as nmp
from netCDF4 import Dataset
import string
from os.path import exists,basename

iconv=2

if   iconv == 1:
    cv_lon = 'lon'
    cv_lat = 'lat'

elif iconv == 2:
    # ERA5...
    cv_lon = 'longitude'
    cv_lat = 'latitude'

rmiss = -9999.

def __chck4f__(cf, script_name=''):
    cmesg = 'ERROR: File '+cf+' does not exist !!!'
    if script_name != '': cmesg = 'ERROR in script '+script_name+': File '+cf+' does not exist !!!'
    if not exists(cf):
        print cmesg ; sys.exit(0)
    else:
        print ' *** will open file '+cf


if len(sys.argv) != 4:
    print 'Usage: '+sys.argv[0]+' <IN_FILE.nc> <FIELD_NAME> <LSM_FILE.nc>'
    sys.exit(0)

cf_in  = sys.argv[1]
cv_in  = sys.argv[2]
cf_lsm = sys.argv[3]

#cf_out = basename(string.replace(cf_in, '_'+cv_in+'_', '_'+cv_in+'masked_'))
cf_out = basename(string.replace(cf_in, '.nc', '_masked.nc'))

__chck4f__(cf_in)
__chck4f__(cf_lsm)

print '\n *** Will generate file '+cf_out+' !\n'

#sys.exit(0)

print ' *** Opening land-sea mask file... '+cf_lsm
f_lsm_in = Dataset(cf_lsm)
xmask = f_lsm_in.variables['lsm'][:,:]
f_lsm_in.close()
print ' *** Land-sea mask read!\n'
idx_land = nmp.where(xmask < 0.5)
dim_lsm = xmask.shape
del xmask


# First need time length:
f_in_in = Dataset(cf_in)
vlon     = f_in_in.variables[cv_lon][:]
cunt_lon = f_in_in.variables[cv_lon].units
clnm_lon = f_in_in.variables[cv_lon].long_name
#print 'LONGITUDE: ', cunt_lon, clnm_lon

# Extracting the longitude 1D array:
vlat     = f_in_in.variables[cv_lat][:]
cunt_lat = f_in_in.variables[cv_lat].units
clnm_lat = f_in_in.variables[cv_lat].long_name
#print 'LATITUDE: ', cunt_lat, clnm_lat

# Extracting time 1D array:
vtime     = f_in_in.variables['time'][:]
cunt_time = f_in_in.variables['time'].units
ccal_time = f_in_in.variables['time'].calendar
#print 'TIME: ', cunt_time, '\n'
f_in_in.close()



Nt = len(vtime)

print 'Nt = ', Nt

for jt in range(Nt):

    print ' *** jt = ', jt


    if jt == 0:
        f_in_in = Dataset(cf_in)
        cunt_in = f_in_in.variables[cv_in].units
        clnm_in = f_in_in.variables[cv_in].long_name
    # Reading field at time jt:
    xfin     = f_in_in.variables[cv_in][jt,:,:]
    if jt == Nt-1: f_in_in.close()


    # Checking dimensions
    if jt == 0:
        dim_in  = xfin.shape
        if dim_in != dim_lsm:
            print 'Shape problem!!!'; print dim_in , dim_lsm
        print '\n'
        (nj,ni) = dim_in
        print 'ni, nj, nt = ', ni, nj, Nt
        xfout = nmp.zeros((nj,ni))


    # Masking land on field
    # ~~~~~~~~~~~~~~~~~~~~~
    xfout[:,:]      = xfin[:,:]
    xfout[idx_land] = rmiss


    # Creating output file
    # ~~~~~~~~~~~~~~~~~~~~
    if jt == 0:
        f_out = Dataset(cf_out, 'w', format='NETCDF4')

        # Dimensions:
        f_out.createDimension(cv_lon, ni)
        f_out.createDimension(cv_lat, nj)
        f_out.createDimension('time', None)

        # Variables
        id_lon = f_out.createVariable(cv_lon, 'f4',(cv_lon,),               zlib=True)
        id_lat = f_out.createVariable(cv_lat, 'f4',(cv_lat,),               zlib=True)
        id_tim = f_out.createVariable('time', 'f4',('time',),               zlib=True)
        id_out  = f_out.createVariable(cv_in, 'f4',('time',cv_lat,cv_lon,), zlib=True, fill_value=rmiss)

        # Attributes
        id_tim.units    = cunt_time
        id_tim.calendar = ccal_time

        id_lat.units         = cunt_lat
        id_lat.long_name     = clnm_lat
        #id_lat.standard_name = csnm_lat

        id_lon.units         = cunt_lon
        id_lon.long_name     = clnm_lon
        #id_lon.standard_name = csnm_lon

        id_out.units     = cunt_in
        id_out.long_name = clnm_in
        #id_out.code  = '133'
        #id_out.table = '128'

        f_out.About = cv_in+' masked and deflated with "mask_field_ecmwf.py" of AeroBulk. [https://github.com/brodeau/aerobulk].'

        # Filling variables:
        id_lat[:] = vlat[:]
        id_lon[:] = vlon[:]

    id_tim[jt]     = vtime[jt]
    id_out[jt,:,:] = xfout[:,:]

    if jt == Nt-1: f_out.close()

print 'Bye!'
