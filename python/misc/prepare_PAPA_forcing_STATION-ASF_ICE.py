#!/usr/bin/env python
#
# L. Brodeau, August 2019

import sys
import numpy as nmp
from netCDF4 import Dataset
import string
from os.path import exists,basename

l_treat_miss_val = False

#cname = 'PAPA' ; ccoor = '50N-145W' ; clon_in = 'lon' ; clat_in = 'lat'
cname = 'ERA5' ; ccoor = '85p5N_54W' ; clon_in = 'longitude' ; clat_in = 'latitude'

#cdir_in = '/home/laurent/PAPA_2012-2018'
cdir_in = '/MEDIA/data/STATION_ASF/ICE'

#cf_in = cdir_in+'/ERA5_station_arctic_1h_2018_1x1.nc'

fext='.nc'

#list_vari = [  'WU_422'     ,  'WV_423'     ,  'WS_401'     ,    'AT_21'       ,   'RH_910'    ,    'Ql_136'   ,    'RD_495'    ,    'BP_915'   ,    'T_25'      ,     'S_41'     , 'U_320'        , 'V_321'         ] #, 'RN_485'         ]
#list_varo = [  'u_air'      ,  'v_air'      ,  'wndspd'     ,    't_air'       ,   'rh_air'    ,    'rad_lw'   ,    'rad_sw'    ,    'slp'      ,     'sst'      ,     'sss'      , 'ssu'          , 'ssv'           ] #, 'rain'           ]
#list_real = [  'atm'        ,  'atm'        ,  'atm'        ,    'atm'         ,   'atm'       ,    'atm'        ,  'atm'        ,   'atm'      ,      'oce'     ,     'oce'      , 'oce'          , 'oce'           ] #, 'atm'            ]
#list_rmlt = [    1.         ,    1.         ,    1.         ,      1.          ,     1.        ,      1.         ,    1.         ,     100.     ,        1.      ,       1.       ,  0.01          ,   0.01          ] #,   1.             ]
#list_rofs = [    0.         ,    0.         ,    0.         ,    273.15        ,     0.        ,      0.         ,    0.         ,     0.       ,        0.      ,       0.       ,  0.            ,   0.            ] #,   1.             ]

list_vari = [      'istl1'     ] 
list_varo = [       'sit'      ] 
list_real = [       'ice'      ] 
list_rmlt = [         1.       ] 
list_rofs = [         0.       ] 

cv_lon = 'nav_lon'
cv_lat = 'nav_lat'
cv_tim = 'time_counter'


def __chck4f__(cf, script_name=''):
    cmesg = 'ERROR: File '+cf+' does not exist !!!'
    if script_name != '': cmesg = 'ERROR in script '+script_name+': File '+cf+' does not exist !!!'
    if not exists(cf):
        print cmesg ; sys.exit(0)
    else:
        print ' *** will open file '+cf




if len(sys.argv) != 3:
    print 'Usage: '+sys.argv[0]+' <file_in> <year>\n'
    sys.exit(0)
cf_in = sys.argv[1]
cyear = sys.argv[2]

#cf_atm = 'Station_'+cname+'_'+ccoor+'_atm_hourly_y'+cyear+'.nc'
cf_ice = 'Station_'+cname+'_'+ccoor+'_ice_hourly_y'+cyear+'.nc'


# First checking if all files are here:
#for cf in list_file: __chck4f__(cdir_in+'/'+cf+fext)
#print '\n ** All files are here, good!\n'


iv = 0
for cv_in in list_vari:

    #cf_in = cdir_in+'/'+list_file[iv]+fext

    rmlt = list_rmlt[iv]
    rofs = list_rofs[iv]
    
    print '\n\n====================================================='
    print '\n *** Doing variable "'+cv_in+'" in file "'+cf_in+'" !'

    #####################################################
    f_in = Dataset(cf_in)

    if iv == 0:
        # Extracting the longitude 1D array:
        vlon     = f_in.variables[clon_in][:]
        cunt_lon = f_in.variables[clon_in].units
        # Extracting the latitude 1D array:
        vlat     = f_in.variables[clat_in][:]
        cunt_lat = f_in.variables[clat_in].units
        # Extracting time 1D array:
        vtime     = f_in.variables['time'][:]
        cunt_time = f_in.variables['time'].units
        Nt0       = len(vtime)

    # How does the variable looks like:
    #xfin0    = f_in.variables[cv_in][:,:,:,:] ; #LOLO: there is this 'depth'
    xfin0    = f_in.variables[cv_in][:,:,:] ; #LOLO: there is this 'depth'
    cunt_in  = f_in.variables[cv_in].units
    clnm_in  = f_in.variables[cv_in].long_name

    (Nt, nj,ni) = nmp.shape(xfin0)
    if Nt != Nt0:             print 'ERROR #1 / (Nt != Nt0)', Nt, Nt0 ;       sys.exit(1)
    if (nj,ni) != (1,1): print 'ERROR #2 / ((nj,ni) != (1,1))'; sys.exit(2)
    print 'ni, nj, Nt = ', ni, nj, Nt
    vfin = nmp.zeros(Nt)
    vfin[:] = xfin0[:,0,0]
    f_in.close()
    del xfin0
    #####################################################



    if l_treat_miss_val:
        ########################################################
        # Finding missing values and doing linear interpolation
        #########################################################
        # First, treat single missing measures:
        for jt in range(1,Nt-1):
            if vfin[jt] > 10000. and vfin[jt-1] < 10000. and vfin[jt+1] < 10000.: vfin[jt] = 0.5*(vfin[jt-1] + vfin[jt+1])
        #
        ( idx_miss, ) = nmp.where(vfin > 10000.)
        nmiss = len(idx_miss)
        if nmiss > 0:
            print ' idx_miss =', idx_miss
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


    # FOR C1D of NEMO, on which STATION-ASF is built upon, time series are 3x3 in space...
    ni = 3
    nj = 3
    
    if iv == 0:

        # Creating output file for ocean:
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        f_ice = Dataset(cf_ice, 'w', format='NETCDF4')
        # Dimensions:
        f_ice.createDimension('x'   , ni  )
        f_ice.createDimension('y'   , nj  )
        f_ice.createDimension(cv_tim, None)
        # Variables
        ido_lon = f_ice.createVariable(cv_lon, 'f4', ('y','x',), zlib=True)
        ido_lat = f_ice.createVariable(cv_lat, 'f4', ('y','x',), zlib=True)
        ido_tim = f_ice.createVariable(cv_tim, 'f4', (cv_tim,) , zlib=True)
        # Attributes
        ido_tim.units     = cunt_time
        ido_lat.units     = cunt_lat
        ido_lon.units     = cunt_lon
        f_ice.About = 'prepare_PAPA_forcing_STATION-ASF_ICE.py of AeroBulk for NEMO/STATION-ASF test-case. Gaps in time-series are filled by means of linear interpolation.'        
        # Filling variables:
        ido_lat[:,:] = vlat[0]
        ido_lon[:,:] = vlon[0]
        for jt in range(Nt): ido_tim[jt] = vtime[jt]

    if list_real[iv] == 'atm':
        id_atm       = f_forc.createVariable(list_varo[iv], 'f4',(cv_tim,'y','x',), zlib=True) ;#, fill_value=rmiss)
        id_atm.units = cunt_in    
        if cv_in in ['AT_21']:  id_atm.units = 'K'
        if cv_in in ['BP_915']: id_atm.units = 'Pa'
        id_atm.long_name = clnm_in            
        for jt in range(Nt): id_atm[jt,:,:] = rmlt*vfin[jt] + rofs
            
    elif list_real[iv] == 'oce':
        id_ice       = f_ice.createVariable(list_varo[iv], 'f4',(cv_tim,'y','x',), zlib=True) ;#, fill_value=rmiss)
        id_ice.units = cunt_in
        if cv_in in ['U_320','V_321']: id_ice.units = 'm/s'
        id_ice.long_name = clnm_in        
        for jt in range(Nt): id_ice[jt,:,:] = rmlt*vfin[jt] + rofs

    elif list_real[iv] == 'ice':
        id_ice       = f_ice.createVariable(list_varo[iv], 'f4',(cv_tim,'y','x',), zlib=True) ;#, fill_value=rmiss)
        id_ice.units = cunt_in
        if cv_in in ['U_320','V_321']: id_ice.units = 'm/s'
        id_ice.long_name = clnm_in        
        for jt in range(Nt): id_ice[jt,:,:] = rmlt*vfin[jt] + rofs


    iv = iv + 1


# Missing fields for ATM:
#for cv in ['rain','snow']:
#    print '\n Creating empty variable '+cv+' into '+cf_atm+' ...'
#    id_atm  = f_forc.createVariable(cv, 'f4',(cv_tim,'y','x',), zlib=True)
#    id_atm.units = '...'
#    id_atm.long_name = 'MISSING!'
#    for jt in range(Nt): id_atm[jt,:,:] = 0.
    
f_ice.close()


## Missing fields for OCE:
#for cv in ['ssh','e3t_m','frq_m']:
#    print '\n Creating empty variable '+cv+' into '+cf_oce+' ...'
#    id_oce  = f_oce.createVariable(cv, 'f4',(cv_tim,'y','x',), zlib=True)
#    id_oce.units = '...'
#    id_oce.long_name = 'MISSING!'
#    for jt in range(Nt): id_oce[jt,:,:] = 0.

#f_oce.close()



del vlon, vlat, vtime



sys.exit(0)


# Daily PRECIP:

# Reading
cf_in = cdir_in+'/'+'rain50n145w_dy.cdf' ; cv_in = 'RN_485' ; cv_out = 'precip' ; rmlt=1./3600. ; rofs=0.
__chck4f__(cf_in)
f_in = Dataset(cf_in)
# Extracting the longitude 1D array:
vlon     = f_in.variables['lon'][:]
cunt_lon = f_in.variables['lon'].units
# Extracting the latitude 1D array:
vlat     = f_in.variables['lat'][:]
cunt_lat = f_in.variables['lat'].units
# Extracting time 1D array:
vtime     = f_in.variables['time'][:]
cunt_time = f_in.variables['time'].units
Nt0       = len(vtime)
# How does the variable looks like:
xfin0    = f_in.variables[cv_in][:,:,:,:] ; #LOLO: there is this 'depth'
#cunt_in  = f_in.variables[cv_in].units
clnm_in  = f_in.variables[cv_in].long_name
(Nt, nk,nj,ni) = nmp.shape(xfin0)
if Nt != Nt0:             print 'ERROR #1 / (Nt != Nt0)', Nt, Nt0 ;       sys.exit(1)
if (nj,ni) != (1,1) or nk > 2: print 'ERROR #2 / ((nk,nj,ni) != (1,1,1))'; sys.exit(2)
print 'ni, nj, nk, Nt = ', ni, nj, nk, Nt
vfin = nmp.zeros(Nt)
vfin[:] = xfin0[:,0,0,0]
f_in.close()
# FOR C1D of NEMO, on which STATION-ASF is built upon, time series are 3x3 in space...
ni = 3
nj = 3






# Writing
cf_atm = 'Station_'+cname+'_'+ccoor+'_precip_daily_y'+cyear+'.nc'
f_atm = Dataset(cf_atm, 'w', format='NETCDF4')
# Dimensions:
f_atm.createDimension('x'   , ni  )
f_atm.createDimension('y'   , nj  )
f_atm.createDimension(cv_tim, None)
# Variables
ida_lon = f_atm.createVariable(cv_lon, 'f4', ('y','x',), zlib=True)
ida_lat = f_atm.createVariable(cv_lat, 'f4', ('y','x',), zlib=True)
ida_tim = f_atm.createVariable(cv_tim, 'f4', (cv_tim,) , zlib=True)
# Attributes
ida_tim.units     = cunt_time
ida_lat.units     = cunt_lat
ida_lon.units     = cunt_lon
f_atm.About = 'Created by L. Brodeau for NEMO/STATION-ASF test-case. Gaps in time-series are filled by means of linear interpolation.'        
# Filling variables:
ida_lat[:,:] = vlat[0]
ida_lon[:,:] = vlon[0]
for jt in range(Nt): ida_tim[jt] = vtime[jt]
#
id_atm       = f_atm.createVariable(cv_out, 'f4',(cv_tim,'y','x',), zlib=True) ;#, fill_value=rmiss)
id_atm.units = 'kg/m^2/s'    
id_atm.long_name = clnm_in            
for jt in range(Nt): id_atm[jt,:,:] = rmlt*vfin[jt] + rofs
#
id_atm       = f_atm.createVariable('snow', 'f4',(cv_tim,'y','x',), zlib=True) ;#, fill_value=rmiss)
id_atm.units = 'kg/m^2/s'
id_atm.long_name = 'Solid precipitation'
for jt in range(Nt): id_atm[jt,:,:] = 0.

f_atm.close()

print 'Bye!'
