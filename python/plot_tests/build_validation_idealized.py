#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-

# Get the "idealized output" of "do_buoy_test_series.sh", run with "IDEALIZED" forcing
# to generate the validation reference for each flux component: a mean value, and a lower and upper value not to go beyond...

import sys
from os import path as path
import math
import numpy as nmp
from netCDF4 import Dataset,num2date
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

cf_out = 'VALIDATION_IDEALIZED.nc'

# Algorithms to include in the diagnostic
l_alg = [ 'andreas' , 'coare3p0' ,  'coare3p6',  'ecmwf',  'ncar' ]
l_clr = [  '#ffed00', 'orange'   , '#008ab8'  ,     '0.4'  , 'pink'  ] ; # colors to differentiate algos on the plot
l_wdt = [     3     ,      2     ,      2     ,       1    ,  2   ] ; # line-width to differentiate algos on the plot
l_lst = [    '-'    ,    '--'    ,     '-'    ,      '--'  , '-'  ] ; # line-style
nba     =    len(l_alg)


# Variables to work with:
l_var = [   'Qlat'   ,    'Qsen'   ,    'Qlw'   ,   'Tau'   ]


dir_figs='.'
size_fig=(13,8)
fig_ext='png'

#clr_red = '#AD0000'
#clr_sat = '#ffed00'
#clr_mod = '#008ab8'

rDPI=100.


nb_algos = len(l_alg)


# Getting arguments:
l_forc_f = False
narg     = len(sys.argv)
if not narg in [2,3]:
    print('Usage: '+sys.argv[0]+' <DIR_OUT_IDEALIZED> (<file_USED_for_idealized_forcing>'); sys.exit(0)
cdir_rslt = sys.argv[1]
if narg==3:
    l_forc_f = True
    cfil_forc = sys.argv[2]
    cf_out = str.replace(cf_out, '.nc', '+forcing.nc')


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Populating and checking existence of files to be read
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def chck4f(cf):
    cmesg = 'ERROR: File '+cf+' does not exist !!!'
    if not path.exists(cf): print(cmesg) ; sys.exit(0)



if l_forc_f: chck4f(cfil_forc)

cf_in = []
for ja in range(nb_algos):
    cfi = cdir_rslt+'/IDEALIZED_'+l_alg[ja]+'.nc'
    chck4f(cfi)
    cf_in.append(cfi)
print('Files we are goin to use:')
for ja in range(nb_algos): print(cf_in[ja])
#-----------------------------------------------------------------


# Getting time array from the first file:
id_in = Dataset(cf_in[0])
vt = id_in.variables['time'][:]
cunit_t = id_in.variables['time'].units ; print(' "time_counter" is in "'+cunit_t+'"')
id_in.close()
Nt = len(vt)

vtime = nmp.zeros(Nt); vtime[:]  = vt[:]
#vtime = num2date(vt, units=cunit_t) ; # something understandable!
#vtime = vtime.astype(dtype='datetime64[D]')


print('\n *** Number of air-sea conditions tested => Nt = '+str(Nt)+'\n')



nb_var = len(l_var)

xF  = nmp.zeros((Nt,nb_algos,nb_var))
xMn = nmp.zeros((Nt,nb_var))    ; # Mean value / algorithms
xSD = nmp.zeros((Nt,nb_var))    ; # Standard Deviation between algorithms


for ja in range(nb_algos):
    print('\n *** Algo = '+l_alg[ja])
    id_in = Dataset(cf_in[ja])
    
    for jv in range(nb_var):
        print('    => reading variable: '+l_var[jv]+' !')
        #
        xF[:,ja,jv] = id_in.variables[l_var[jv]][:] # only the center point of the 3x3 spatial domain!

    id_in.close()


for jv in range(nb_var):
    vm = nmp.mean(xF[:,:,jv], axis=1)
    xMn[:,jv] = vm
    for jt in range(Nt):
        rr         = xF[jt,:,jv]-vm[jt]
        xSD[jt,jv] = nmp.sqrt( nmp.mean( rr*rr ) )


rtol = 2.2  ; # tolerance in terms of number of times the SD !


# Testing if the tolerance of "rtol" standard deviation is good:
for jv in range(nb_var):
    for ja in range(nb_algos):
        ll_a = nmp.any( xF[:,ja,jv] > xMn[:,jv]+rtol*xSD[:,jv]  )
        ll_b = nmp.any( xF[:,ja,jv] < xMn[:,jv]-rtol*xSD[:,jv]  )


if ll_a: print('Well, seems like the choice of '+str(rtol)+' standard deviation is no sufficient ABOVE!')
if ll_b: print('Well, seems like the choice of '+str(rtol)+' standard deviation is no sufficient BELOW!')
if ll_a or ll_b: sys.exit(0)

ii=Nt/300
ib=max(ii-ii%10,1)
xticks_d=int(30*ib)

rat = 100./float(rDPI)
params = { 'font.family':'Open Sans',
           'font.size':       int(15.*rat),
           'legend.fontsize': int(15.*rat),
           'xtick.labelsize': int(15.*rat),
           'ytick.labelsize': int(15.*rat),
           'axes.labelsize':  int(16.*rat)
}
mpl.rcParams.update(params)
font_inf = { 'fontname':'Open Sans', 'fontweight':'normal', 'fontsize':18.*rat }
font_x   = { 'fontname':'Open Sans', 'fontweight':'normal', 'fontsize':15.*rat }


# Now we compare output variables from bulk algorithms between them:

for jv in range(nb_var):
        
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fig = plt.figure(num = jv, figsize=size_fig, facecolor='w', edgecolor='k')
    ax1 = plt.axes([0.08, 0.25, 0.9, 0.7])
    #ax1.set_xticks(vtime[::xticks_d])
    #ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
    #plt.xticks(rotation='60', **font_x)

    for ja in range(nb_algos):
        plt.plot(vtime, xF[:,ja,jv], color=l_clr[ja], linestyle=l_lst[ja], linewidth=l_wdt[ja], label=l_alg[ja], zorder=10+ja)


    plt.plot(vtime, xMn[:,jv], color='k', linestyle='-.', linewidth=2, label='Mean', zorder=100)
    plt.plot(vtime, xMn[:,jv]+rtol*xSD[:,jv], color='r', linestyle='-.', linewidth=2, label='Mean + '+str(rtol)+'*SD', zorder=100)
    plt.plot(vtime, xMn[:,jv]-rtol*xSD[:,jv], color='b', linestyle='-.', linewidth=2, label='Mean - '+str(rtol)+'*SD', zorder=100)
    

    #ax1.set_ylim(L_VMIN[jv], L_VMAX[jv]) ; ax1.set_xlim(vtime[0],vtime[Nt-1])
    #plt.ylabel(L_VARL[jv]+' ['+L_VUNT[jv]+']')

    ax1.grid(color='k', linestyle='-', linewidth=0.3)
    plt.legend(loc='best', ncol=1, shadow=True, fancybox=True)

    plt.savefig(l_var[jv]+'.'+fig_ext, dpi=int(rDPI), transparent=False)
    plt.close(jv)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


if l_forc_f: idf_f = Dataset(cfil_forc) ; # forcing that was used to generate the results
idf_i = Dataset(cf_in[0])  ; # One of the input file of the results
idf_o = Dataset(cf_out, 'w')
# cfil_forc





if l_forc_f:

    print('\n *** Wild add all forcing fields (ocean+atmo), that were used to generate the results, into '+cf_out+' !\n')
    
    # copy global attributes all at once via dictionary
    idf_o.setncatts(idf_f.__dict__)

    # copy dimensions
    for name, dimension in idf_f.dimensions.items():
        idf_o.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))

    # copy all file data except for the excluded
    for name, variable in idf_f.variables.items():
        x = idf_o.createVariable(name, variable.datatype, variable.dimensions)
        idf_o[name][:] = idf_f[name][:]
        # copy variable attributes all at once via dictionary
        idf_o[name].setncatts(idf_f[name].__dict__)

else:
    print('\n *** NOT adding all forcing fields (ocean+atmo) into '+cf_out+' !\n')
    # copy global attributes all at once via dictionary:
    idf_o.setncatts(idf_i.__dict__)
    # copy dimensions:
    for name, dimension in idf_i.dimensions.items():
        idf_o.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))


for name, variable in idf_i.variables.items():
    #print('name, variable =', name, variable)

    if not l_forc_f:
        if name == 'time':
            x = idf_o.createVariable(name, variable.datatype, variable.dimensions)
            idf_o[name].setncatts(idf_i[name].__dict__)
            x[:] = vtime[:]

    if name in l_var:
        
        jv = l_var.index(name)
        print('\n YES! =>', name, jv, '\n')

        x1 = idf_o.createVariable(name+'_mean', variable.datatype, variable.dimensions)
        idf_o[name+'_mean'].setncatts(idf_i[name].__dict__)
        
        x2 = idf_o.createVariable(name+'_tol', variable.datatype, variable.dimensions)
        idf_o[name+'_tol'].setncatts(idf_i[name].__dict__)


        
        x1[:] = xMn[:,jv]
        x2[:] = xSD[:,jv]



    
    
idf_o.close()
idf_i.close()
if l_forc_f: idf_f.close()
sys.exit(0)





# Saving the validation reference result:

idf_i = Dataset(cf_in[0]) 
idf_o = Dataset('VALIDATION_IDEALIZED.nc', 'w')

# copy global attributes all at once via dictionary
idf_o.setncatts(idf_i.__dict__)

# copy dimensions
for name, dimension in idf_i.dimensions.items():
    idf_o.createDimension(
        name, (len(dimension) if not dimension.isunlimited() else None))
    
for name, variable in idf_i.variables.items():
    #print('name, variable =', name, variable)

    if name == 'time':
        x = idf_o.createVariable(name, variable.datatype, variable.dimensions)
        idf_o[name].setncatts(idf_i[name].__dict__)
        x[:] = vtime[:]
        
    if name in l_var:
        
        jv = l_var.index(name)
        print('\n YES! =>', name, jv, '\n')

        x1 = idf_o.createVariable(name+'_mean', variable.datatype, variable.dimensions)
        idf_o[name+'_mean'].setncatts(idf_i[name].__dict__)
        
        x2 = idf_o.createVariable(name+'_tol', variable.datatype, variable.dimensions)
        idf_o[name+'_tol'].setncatts(idf_i[name].__dict__)


        
        x1[:] = xMn[:,jv]
        x2[:] = xSD[:,jv]

idf_o.close()        
idf_i.close()

