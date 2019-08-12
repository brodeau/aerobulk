#!/usr/bin/env python
#
# L. Brodeau, August 2019

import sys
import numpy as nmp
from netCDF4 import Dataset
import string
from os.path import exists,basename

cdir_in = '/home/laurent/PAPA_2012-2018'

fext='.cdf'

list_file = [ 'w50n145w_hr' , 'airt50n145w_hr' , 'rh50n145w_hr', 'lw50n145w_hr', 'rad50n145w_hr', 'bp50n145w_hr', 'sst50n145w_hr' ]
list_vari = [  'WS_401'     ,    'AT_21'       ,   'RH_910'    ,    'Ql_136'   ,    'RD_495'    ,    'BP_915'   ,    'T_25'       ]
list_varo = [  'wndspd'     ,    't_air'       ,   'rh_air'    ,    'rad_lw'   ,    'rad_sw'    ,      'slp'    ,     'sst'       ]

cf_out = 'Station_PAPA_50N-145W_2012-2018.nc4'

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

    cf_in = cdir_in+'/'+cf+fext
    cv_in = list_vari[iv]

    print '\n\n====================================================='
    __chck4f__(cf_in)
    print '\n *** Doing variable "'+cv_in+'" in file "'+cf_in+'" !'

    #####################################################
    f_in = Dataset(cf_in)

    if iv == 0:
        # Extracting the longitude 1D array:
        vlon     = f_in.variables[cv_lon][:]
        cunt_lon = f_in.variables[cv_lon].units
        # Extracting the latitude 1D array:
        vlat     = f_in.variables[cv_lat][:]
        cunt_lat = f_in.variables[cv_lat].units
        # Extracting time 1D array:
        vtime     = f_in.variables['time'][:]
        cunt_time = f_in.variables['time'].units
        Nt0       = len(vtime)

    # How does the variable looks like:
    xfin0    = f_in.variables[cv_in][:,:,:,:] ; #LOLO: there is this 'depth'
    cunt_in  = f_in.variables[cv_in].units
    clnm_in  = f_in.variables[cv_in].long_name

    (Nt, nk,nj,ni) = nmp.shape(xfin0)
    if Nt != Nt0:             print 'ERROR #1 / (Nt != Nt0)';             sys.exit(1)
    if (nk,nj,ni) != (1,1,1): print 'ERROR #2 / ((nk,nj,ni) != (1,1,1))'; sys.exit(2)
    print 'ni, nj, nk, Nt = ', ni, nj, nk, Nt
    vfin = nmp.zeros(Nt)
    vfin[:] = xfin0[:,0,0,0]
    f_in.close()
    del xfin0
    #####################################################

    
    ########################################################
    # Finding missing values and doing linear interpolation
    #########################################################
    # First, treat single missing measures:
    for jt in range(1,Nt-1):
        if vfin[jt] > 10000. and vfin[jt-1] < 10000. and vfin[jt+1] < 10000.: vfin[jt] = 0.5*(vfin[jt-1] + vfin[jt+1])
    #
    ( idx_miss, ) = nmp.where(vfin > 10000.)
    print ' idx_miss =', idx_miss
    nmiss = len(idx_miss)
    jt0 = idx_miss[0]
    jtN = idx_miss[nmiss-1]
    print ' jt0, jtN = ', jt0, jtN
    ibndmiss = nmp.zeros(nmiss)
    ii = 0
    ibndmiss[ii] = 1
    while ii < nmiss-1:
        ii = ii + 1
        if idx_miss[ii] > idx_miss[ii-1] + 1 :
            ibndmiss[ii-1] = 2
            ibndmiss[ii] = 1
    ibndmiss[nmiss-1] = 2
    print ' ibndmiss =', ibndmiss
    (idx1,) = nmp.where(ibndmiss==1)
    (idx2,) = nmp.where(ibndmiss==2)    
    if len(idx1) != len(idx2): print 'PROBLEM with idx1 and idx2 !'; sys.exit(3)
    # => if only 1 missing point 
    nseg = len(idx1)
    print ' idx1, idx2 =', idx1, idx2
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

    if iv == 0:
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

        # Attributes
        id_tim.units     = cunt_time
        id_lat.units     = cunt_lat
        id_lon.units     = cunt_lon

        f_out.About = 'Created by L. Brodeau for AeroBulk. Gaps in time-series are filled by means of linear interpolation.'
        
        # Filling variables:
        id_lat[:] = vlat[:]
        id_lon[:] = vlon[:]
        for jt in range(Nt): id_tim[jt] = vtime[jt]

        
    id_out  = f_out.createVariable(list_varo[iv], 'f4',('time',cv_lat,cv_lon,), zlib=True) ;#, fill_value=rmiss)
    id_out.units     = cunt_in
    id_out.long_name = clnm_in    
    for jt in range(Nt): id_out[jt,0,0] = vfin[jt]



    iv = iv + 1
        
f_out.close()


print 'Bye!'
