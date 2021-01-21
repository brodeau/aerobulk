#!/usr/bin/env python3
#
# L. Brodeau, June 2019

import sys
import numpy as nmp
from netCDF4 import Dataset
from os.path import exists,basename

import aerobulk_physf as abp

rmiss = -9999.

def __chck4f__(cf, script_name=''):
    cmesg = 'ERROR: File '+cf+' does not exist !!!'
    if script_name != '': cmesg = 'ERROR in script '+script_name+': File '+cf+' does not exist !!!'
    if not exists(cf):
        print(cmesg) ; sys.exit(0)
    else:
        print(' *** will open file '+cf)

if not len(sys.argv) in [2,3]:
    print('Usage: '+sys.argv[0]+' <IN_FILE_D2.nc> (<lsm_file.nc>)')
    sys.exit(0)


cf_d2   = sys.argv[1]

l_mask = False
if len(sys.argv) == 3:
    l_mask = True
    cf_lsm = sys.argv[2]

if l_mask:
    print(' *** Opening land-sea mask file... '+cf_lsm)
    f_lsm_in = Dataset(cf_lsm)
    xmask = f_lsm_in.variables['lsm'][:,:]
    f_lsm_in.close()
    print(' *** Land-sea mask read!\n')
    idx_land = nmp.where(xmask < 0.5)
    del xmask




# First need time length:
f_d2_in = Dataset(cf_d2)
#
list_var = f_d2_in.variables.keys()

# Name of longitude and latitude ?
cv_lon='0'
for cv_lon in [ 'lon', 'longitude', 'LON', 'LONGITUDE', 'nav_lon', 'glamt' ]:
    if cv_lon in list_var: break
if cv_lon=='0': print("PROBLEM: could not find name of longitude variable!"); sys.exit(0)
#
cv_lat='0'
for cv_lat in [ 'lat', 'latitude', 'LAT', 'LATITUDE', 'nav_lat', 'gphit' ]:
    if cv_lat in list_var: break
if cv_lat=='0': print("PROBLEM: could not find name of latitude variable!"); sys.exit(0)
#
cv_time='0'
for cv_time in [ 'time', 'time_counter', 'TIME' ]:
    if cv_time in list_var: break
if cv_time=='0': print("PROBLEM: could not find name of time variable!"); sys.exit(0)
#
cv_d2='0'
for cv_d2 in [ 'd2m', 'D2M' ]:
    if cv_d2 in list_var: break
if cv_d2=='0': print("PROBLEM: could not find name of d2 variable!"); sys.exit(0)
if   cv_d2=='D2M':
    cv_p0='MSL'
    cv_q2='Q2M'
elif cv_d2=='d2m':
    cv_p0='msl'
    cv_q2='q2m'

cunt_lon='unknown' ; clnm_lon='Longitude'
vlon    = f_d2_in.variables[cv_lon][:]
lst_att = f_d2_in.variables[cv_lon].ncattrs()
if 'units'     in lst_att: cunt_lon = f_d2_in.variables[cv_lon].units
if 'long_name' in lst_att: clnm_lon = f_d2_in.variables[cv_lon].long_name
print('LONGITUDE: ', cunt_lon, clnm_lon)
#
cunt_lat='unknown' ; clnm_lat='Latitude'
vlat    = f_d2_in.variables[cv_lat][:]
lst_att = f_d2_in.variables[cv_lat].ncattrs()
if 'units'     in lst_att: cunt_lat = f_d2_in.variables[cv_lat].units
if 'long_name' in lst_att: clnm_lat = f_d2_in.variables[cv_lat].long_name
print('LATITUDE: ', cunt_lat, clnm_lat)
#
cunt_time='unknown'
vtime   = f_d2_in.variables[cv_time][:]
lst_att = f_d2_in.variables[cv_time].ncattrs()
if 'units' in lst_att: cunt_time = f_d2_in.variables[cv_time].units
print('TIME: ', cunt_time, '\n')
#
f_d2_in.close()



cf_p0 = str.replace(cf_d2, cv_d2, cv_p0)
__chck4f__(cf_p0)

cf_q2 = basename(str.replace(cf_d2, cv_d2, cv_q2))
__chck4f__(cf_d2)
if l_mask:
    __chck4f__(cf_lsm)
    cf_q2 = basename(str.replace(cf_q2, '.nc', '_masked.nc'))
print('\n *** Will generate file '+cf_q2+' !\n')




Nt = len(vtime)
#Nt=4
print('Nt = ', Nt)

for jt in range(Nt):

    print(' *** jt = ', jt)
    
        
    # D2M
    # ~~~
    if jt == 0:
        f_d2_in = Dataset(cf_d2)
        cunt_d2 = 'unknown'
        lst_att = f_d2_in.variables[cv_d2].ncattrs()
        if 'units' in lst_att: cunt_d2 = f_d2_in.variables[cv_d2].units
    xd2     = f_d2_in.variables[cv_d2][jt,:,:]
    if jt == Nt-1: f_d2_in.close()
    
    
    # MSL
    # ~~~
    if jt == 0:
        f_p0_in = Dataset(cf_p0)
        cunt_p0 = 'unknown'
        lst_att = f_p0_in.variables[cv_p0].ncattrs()
        if 'units' in lst_att: cunt_p0 = f_p0_in.variables[cv_p0].units
    xp0     = f_p0_in.variables[cv_p0][jt,:,:]
    if jt == Nt-1: f_p0_in.close()
                
    
    # Checking dimensions
    # ~~~~~~~~~~~~~~~~~~~
    if jt == 0:
        dim_d2 = xd2.shape ; dim_p0 = xp0.shape
        if dim_d2 != dim_p0:
            print('Shape problem!!!'); print(dim_d2 , dim_p0)
        print('\n')
        [ nj, ni ] = dim_d2
        print('ni, nj, nt = ', ni, nj, Nt)
        xq2 = nmp.zeros(nj*ni) ; xq2.shape = dim_d2
    
    
    # Building q2
    # ~~~~~~~~~~~
    xq2 = abp.q_air_dp(xd2, xp0)               

    if l_mask: xq2[idx_land] = rmiss

    
    # Creating output file
    # ~~~~~~~~~~~~~~~~~~~~
    if jt == 0:
        f_out = Dataset(cf_q2, 'w', format='NETCDF4')
    
        # Dimensions:
        f_out.createDimension(cv_lon, ni)
        f_out.createDimension(cv_lat, nj)
        f_out.createDimension(cv_time, None)
    
        # Variables
        id_lon = f_out.createVariable(cv_lon,'f4',(cv_lon,),               zlib=True)
        id_lat = f_out.createVariable(cv_lat,'f4',(cv_lat,),               zlib=True)
        id_tim = f_out.createVariable(cv_time,'f4',(cv_time,),               zlib=True)
        
        if l_mask:
            id_q2  = f_out.createVariable(cv_q2, 'f4',(cv_time,cv_lat,cv_lon,), zlib=True, fill_value=rmiss)
        else:
            id_q2  = f_out.createVariable(cv_q2, 'f4',(cv_time,cv_lat,cv_lon,), zlib=True)

        # Attributes
        id_tim.units    = cunt_time
        #id_tim.calendar = ccal_time
    
        id_lat.units         = cunt_lat
        id_lat.long_name     = clnm_lat
        #id_lat.standard_name = csnm_lat
    
        id_lon.units         = cunt_lon
        id_lon.long_name     = clnm_lon
        #id_lon.standard_name = csnm_lon

        id_q2.units = 'kg/kg'
        id_q2.long_name = 'Surface specific humidity at 2m, built from '+cv_d2+' and '+cv_p0
        id_q2.code  = '133'
        id_q2.table = '128'
    
        f_out.About = 'Created with "build_q2_from_d2_slp.py" of AeroBulk, using '+cv_p0+' and '+cv_d2+'. [https://github.com/brodeau/aerobulk]'
    
        # Filling variables:
        id_lat[:] = vlat[:]
        id_lon[:] = vlon[:]
        
    id_tim[jt]     = vtime[jt]
    id_q2[jt,:,:]  = xq2[:,:] 
    
    if jt == Nt-1: f_out.close()    
        
print('Bye!')
