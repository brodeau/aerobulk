#!/usr/bin/env python
#
# L. Brodeau, August 2019

import sys
import numpy as nmp
from netCDF4 import Dataset
import string
from os.path import exists,basename


cdir_in = '/home/laurent/PAPA_2012-2018'


list_file = [ 'airt50n145w_hr.cdf' ]
list_vari = [       'AT_21'        ]
list_varo = [       't_air'        ]

cf_out = 'test.nc4'

iconv=1

if   iconv == 1:
    cv_lon = 'lon'
    cv_lat = 'lat'

elif iconv == 2:
    # ERA5...
    cv_lon = 'longitude'
    cv_lat = 'latitude'

#rmiss = -9999.

def __chck4f__(cf, script_name=''):
    cmesg = 'ERROR: File '+cf+' does not exist !!!'
    if script_name != '': cmesg = 'ERROR in script '+script_name+': File '+cf+' does not exist !!!'
    if not exists(cf):
        print cmesg ; sys.exit(0)
    else:
        print ' *** will open file '+cf




iv = 0
for cf in list_file:

    cf_in = cdir_in+'/'+cf
    cv_in = list_vari[iv]

    __chck4f__(cf_in)
    print '\n *** Doing variable "'+cv_in+'" in file "'+cf_in+'" !'

    #####################################################
    f_in_in = Dataset(cf_in)

    # Extracting the longitude 1D array:
    vlon     = f_in_in.variables[cv_lon][:]
    cunt_lon = f_in_in.variables[cv_lon].units

    # Extracting the latitude 1D array:
    vlat     = f_in_in.variables[cv_lat][:]
    cunt_lat = f_in_in.variables[cv_lat].units

    # Extracting time 1D array:
    vtime     = f_in_in.variables['time'][:]
    cunt_time = f_in_in.variables['time'].units

    # How does the variable looks like:
    xfin0    = f_in_in.variables[cv_in][:,:,:,:] ; #LOLO: there is this 'depth'
    cunt_in  = f_in_in.variables[cv_in].units
    clnm_in  = f_in_in.variables[cv_in].long_name

    (Nt, nk,nj,ni) = nmp.shape(xfin0)
    if Nt != len(vtime):      print 'ERROR #1'; sys.exit(1)
    if (nk,nj,ni) != (1,1,1): print 'ERROR #2'; sys.exit(2)
    print 'ni, nj, nk, Nt = ', ni, nj, nk, Nt
    vfin = nmp.zeros(Nt)
    vfin[:] = xfin0[:,0,0,0]
    f_in_in.close()
    del xfin0
    #####################################################

    
    ########################################################
    # Finding missing values and doing linear interpolation
    #########################################################
    ( idx_miss, ) = nmp.where(vfin > 1000.)
    #print ' idx_miss =', idx_miss
    nmiss = len(idx_miss)
    jt0 = idx_miss[0]
    jtN = idx_miss[nmiss-1]
    #print ' jt0, jtN = ', jt0, jtN
    ibndmiss = nmp.zeros(nmiss)
    ii = 0
    ibndmiss[ii] = 1
    while ii < nmiss-1:
        ii = ii + 1
        if idx_miss[ii] > idx_miss[ii-1] + 1 :
            ibndmiss[ii-1] = 2
            ibndmiss[ii] = 1
    ibndmiss[nmiss-1] = 2
    #print ' ibndmiss =', ibndmiss
    (idx1,) = nmp.where(ibndmiss==1)
    (idx2,) = nmp.where(ibndmiss==2)
    if len(idx1) != len(idx2): print 'PROBLEM with idx1 and idx2 !'; sys.exit(3)
    nseg = len(idx1)
    #print ' idx1, idx2 =', idx1, idx2
    print '\nThere are '+str(nseg)+' segments of missing values:'
    for isg in range(nseg):
        jt1 = idx_miss[idx1[isg]] ; jta = jt1 - 1
        jt2 = idx_miss[idx2[isg]] ; jtb = jt2 + 1
        print ' *** Seg. #'+str(isg)+' starts at '+str(jt1)+' and stops at '+str(jt2)
        #
        # Linear interpolation to fill the gap defined by the segment:
        print ' => linear interpolation to fill the gap...'
        print ' BEFORE: ', vfin[jta:jtb+1]
        rslp = (vfin[jtb] - vfin[jta]) / (vtime[jtb] - vtime[jta])
        for jt in range(jt1,jt2+1): vfin[jt] = vfin[jta] + (vtime[jt] - vtime[jta])*rslp
        print ' AFTER: ', vfin[jta:jtb+1], '\n'

    ########################################################

    # Creating output file
    # ~~~~~~~~~~~~~~~~~~~~
    f_out = Dataset(cf_out, 'w', format='NETCDF4')

    # Dimensions:
    f_out.createDimension(cv_lon, ni)
    f_out.createDimension(cv_lat, nj)
    f_out.createDimension('time', None)

    # Variables
    id_lon = f_out.createVariable(cv_lon, 'f4',(cv_lon,),               zlib=True)
    id_lat = f_out.createVariable(cv_lat, 'f4',(cv_lat,),               zlib=True)
    id_tim = f_out.createVariable('time', 'f4',('time',),               zlib=True)
    id_out  = f_out.createVariable(list_varo[iv], 'f4',('time',cv_lat,cv_lon,), zlib=True) ;#, fill_value=rmiss)

    # Attributes
    id_tim.units     = cunt_time
    id_lat.units     = cunt_lat
    id_lon.units     = cunt_lon
    id_out.units     = cunt_in
    id_out.long_name = clnm_in    
    #f_out.About = cv_in+' masked and deflated with "mask_field_ecmwf.py" of AeroBulk. [https://github.com/brodeau/aerobulk].'

    # Filling variables:
    id_lat[:] = vlat[:]
    id_lon[:] = vlon[:]
    
    
    for jt in range(Nt):
        id_tim[jt]     = vtime[jt]
        id_out[jt,0,0] = vfin[jt]

        
f_out.close()
print 'Bye!'
