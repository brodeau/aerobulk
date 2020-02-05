#!/usr/bin/env python
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-

# Post-diagnostic of STATION_ASF /  L. Brodeau, 2019

import sys
from os import path as path
#from string import replace
import math
import numpy as nmp
#import scipy.signal as signal
from netCDF4 import Dataset
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
#from string import find
#import warnings
#warnings.filterwarnings("ignore")
#import time

#import barakuda_plot as bp
#import barakuda_tool as bt

reload(sys)
sys.setdefaultencoding('utf8')

cy1     = '2016' ; # First year
cy2     = '2018' ; # Last year

jt0 = 0
jt0 = 17519


dir_figs='.'
size_fig=(13,7)
fig_ext='png'

clr_red = '#AD0000'
clr_sat = '#ffed00'
clr_mod = '#008ab8'

rDPI=200.

L_ALGOS = [ 'COARE3p6' , 'ECMWF'   , 'NCAR' ]
l_xtrns = [ '-noskin'  , '-noskin' ,  ''    ] ; # string to add to algo name (L_ALGOS) to get version without skin params turned on
l_color = [  '#ffed00' , '#008ab8' , '0.4'  ] ; # colors to differentiate algos on the plot
l_width = [     3      ,    2      ,  1     ] ; # line-width to differentiate algos on the plot
l_style = [    '-'     ,   '-'     , '--'   ] ; # line-style

L_VNEM  = [   'qla'     ,     'qsb'     ,     'qt'     ,   'qlw'     ,  'taum'     ,    'dt_skin'         ]
L_VARO  = [   'Qlat'    ,    'Qsen'     ,     'Qnet'   ,   'Qlw'     ,  'Tau'      ,    'dT_skin'         ] ; # name of variable on figure
L_VARL  = [ r'$Q_{lat}$', r'$Q_{sens}$' , r'$Q_{net}$' , r'$Q_{lw}$' , r'$|\tau|$' , r'$\Delta T_{skin}$' ] ; # name of variable in latex mode
L_VUNT  = [ r'$W/m^2$'  , r'$W/m^2$'    , r'$W/m^2$'   , r'$W/m^2$'  , r'$N/m^2$'  ,      'K'             ]
L_VMAX  = [     75.     ,     75.       ,    800.      ,     25.     ,    1.2      ,      -0.7            ]
L_VMIN  = [   -250.     ,   -125.       ,   -400.      ,   -150.     ,    0.       ,       0.7            ]
L_ANOM  = [   True      ,    True       ,    True      ,    True     ,   True      ,      False           ]

#L_VNEM  = [   'qlw'      ]
#L_VARO  = [   'Qlw'      ] ; # name of variable on figure
#L_VARL  = [ r'$Q_{lw}$'  ] ; # name of variable in latex mode
#L_VUNT  = [ r'$W/m^2$'   ]
#L_VMAX  = [     25.      ]
#L_VMIN  = [   -150.      ]
#L_ANOM  = [    True      ]



nb_algos = len(L_ALGOS) ; print(nb_algos)

# Getting arguments:
narg = len(sys.argv)
if narg != 2:
    print 'Usage: '+sys.argv[0]+' <DIR_OUT_SASF>'; sys.exit(0)
cdir_data = sys.argv[1]



# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Populating and checking existence of files to be read
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def chck4f(cf):
    cmesg = 'ERROR: File '+cf+' does not exist !!!'
    if not path.exists(cf): print cmesg ; sys.exit(0)

###cf_in = nmp.empty((), dtype="S10")
cf_in = [] ; cf_in_ns = []
for ja in range(nb_algos):
    cfi = cdir_data+'/output/'+'STATION_ASF-'+L_ALGOS[ja]+'_1h_'+cy1+'0101_'+cy2+'1231_gridT.nc'
    chck4f(cfi)
    cf_in.append(cfi)
# Same but without skin params:
for ja in range(nb_algos):
    cfi = cdir_data+'/output/'+'STATION_ASF-'+L_ALGOS[ja]+l_xtrns[ja]+'_1h_'+cy1+'0101_'+cy2+'1231_gridT.nc'
    chck4f(cfi)
    cf_in_ns.append(cfi)
print('Files we are goin to use:')
for ja in range(nb_algos): print(cf_in[ja])
print('   --- same without cool-skin/warm-layer:')
for ja in range(nb_algos): print(cf_in_ns[ja])
#-----------------------------------------------------------------


# Getting time array from the first file:
id_in = Dataset(cf_in[0])
vt = id_in.variables['time_counter'][jt0:]
cunit_t = id_in.variables['time_counter'].units ; print(' "time_counter" is in "'+cunit_t+'"')
id_in.close()
nbr = len(vt)




vtime = nmp.zeros(nbr)

vt = vt + 1036800. + 30.*60. # BUG!??? don't get why false in epoch to date conversion, and yet ncview gets it right!
for jt in range(nbr): vtime[jt] = mdates.epoch2num(vt[jt])

ii=nbr/300
ib=max(ii-ii%10,1)
xticks_d=int(30*ib)

font_inf = { 'fontname':'Open Sans', 'fontweight':'normal', 'fontsize':14 }

nb_var = len(L_VNEM)

xF  = nmp.zeros((nbr,nb_algos))
xFa = nmp.zeros((nbr,nb_algos))


for ctest in ['skin','noskin']:

    for jv in range(nb_var):
        print('\n *** Treating variable: '+L_VARO[jv]+' ('+ctest+') !')

        for ja in range(nb_algos):
            #
            if ctest == 'skin':   id_in = Dataset(cf_in[ja])
            if ctest == 'noskin': id_in = Dataset(cf_in_ns[ja])
            xF[:,ja] = id_in.variables[L_VNEM[jv]][jt0:,1,1] # only the center point of the 3x3 spatial domain!
            if ja == 0: cvar_lnm = id_in.variables[L_VNEM[jv]].long_name
            id_in.close()

        fig = plt.figure(num = jv, figsize=size_fig, facecolor='w', edgecolor='k')

        ax1 = plt.axes([0.07, 0.22, 0.9, 0.75])

        ax1.set_xticks(vtime[::xticks_d])
        ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
        plt.xticks(rotation='60')

        for ja in range(nb_algos):
            plt.plot(vtime, xF[:,ja], '-', color=l_color[ja], linestyle=l_style[ja], linewidth=l_width[ja], label=L_ALGOS[ja], zorder=10+ja)

        ax1.set_ylim(L_VMIN[jv], L_VMAX[jv]) ; ax1.set_xlim(vtime[0],vtime[nbr-1])
        plt.ylabel(L_VARL[jv]+' ['+L_VUNT[jv]+']')

        ax1.grid(color='k', linestyle='-', linewidth=0.3)
        plt.legend(bbox_to_anchor=(0.45, 0.2), ncol=1, shadow=True, fancybox=True)
        ax1.annotate(cvar_lnm+' ('+ctest+')', xy=(0.3, 0.97), xycoords='axes fraction',  bbox={'facecolor':'w', 'alpha':1., 'pad':10}, zorder=50, **font_inf)
        plt.savefig(L_VARO[jv]+'_'+ctest+'.'+fig_ext, dpi=int(rDPI), transparent=False)
        plt.close(jv)



        if L_ANOM[jv]:

            for ja in range(nb_algos): xFa[:,ja] = xF[:,ja] - nmp.mean(xF,axis=1)

            if nmp.sum(xFa[:,:]) == 0.0:
                print('     Well! Seems that for variable '+L_VARO[jv]+', choice of algo has no impact a all!')
                print('          ==> skipping anomaly plot...')

            else:

                # Want a symetric y-range that makes sense for the anomaly we're looking at:
                rmax = nmp.max(xFa) ; rmin = nmp.min(xFa)
                rmax = max( abs(rmax) , abs(rmin) )
                romagn = math.floor(math.log(rmax, 10)) ; # order of magnitude of the anomaly  we're dealing with
                rmlt = 10.**(int(romagn)) / 2.
                yrng = math.copysign( math.ceil(abs(rmax)/rmlt)*rmlt , rmax)
                #print 'yrng = ', yrng ;  #sys.exit(0)

                fig = plt.figure(num = 10+jv, figsize=size_fig, facecolor='w', edgecolor='k')
                ax1 = plt.axes([0.07, 0.22, 0.9, 0.75])

                ax1.set_xticks(vtime[::xticks_d])
                ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
                plt.xticks(rotation='60')

                for ja in range(nb_algos):
                    plt.plot(vtime, xFa[:,ja], '-', color=l_color[ja], linewidth=l_width[ja], label=L_ALGOS[ja], zorder=10+ja)

                ax1.set_ylim(-yrng,yrng) ; ax1.set_xlim(vtime[0],vtime[nbr-1])
                plt.ylabel(L_VARL[jv]+' ['+L_VUNT[jv]+']')
                ax1.grid(color='k', linestyle='-', linewidth=0.3)
                plt.legend(bbox_to_anchor=(0.45, 0.2), ncol=1, shadow=True, fancybox=True)
                ax1.annotate('Anomaly of '+cvar_lnm+' ('+ctest+')', xy=(0.3, 0.97), xycoords='axes fraction',  bbox={'facecolor':'w', 'alpha':1., 'pad':10}, zorder=50, **font_inf)
                plt.savefig(L_VARO[jv]+'_'+ctest+'_anomaly.'+fig_ext, dpi=int(rDPI), transparent=False)
                plt.close(10+jv)




# Difference skin vs noskin:
xFns = nmp.zeros((nbr,nb_algos))

for jv in range(nb_var-1):
    print('\n *** Treating variable: '+L_VARO[jv]+' ('+ctest+') !')

    for ja in range(nb_algos-1):
        id_in = Dataset(cf_in[ja])
        xF[:,ja]   = id_in.variables[L_VNEM[jv]][jt0:,1,1] # only the center point of the 3x3 spatial domain!
        if ja == 0: cvar_lnm = id_in.variables[L_VNEM[jv]].long_name
        id_in.close()
        #
        id_in = Dataset(cf_in_ns[ja])
        xFns[:,ja] = id_in.variables[L_VNEM[jv]][jt0:,1,1] # only the center point of the 3x3 spatial domain!
        if ja == 0: cvar_lnm = id_in.variables[L_VNEM[jv]].long_name
        id_in.close()

        xFa[:,ja] = xF[:,ja] - xFns[:,ja]  ; # difference!


    # Want a symetric y-range that makes sense for the anomaly we're looking at:
    rmax = nmp.max(xFa) ; rmin = nmp.min(xFa)
    rmax = max( abs(rmax) , abs(rmin) )
    romagn = math.floor(math.log(rmax, 10)) ; # order of magnitude of the anomaly  we're dealing with
    rmlt = 10.**(int(romagn)) / 2.
    yrng = math.copysign( math.ceil(abs(rmax)/rmlt)*rmlt , rmax)
    print 'yrng = ', yrng ;  #sys.exit(0)




    for ja in range(nb_algos-1):

        calgo = L_ALGOS[ja]

        if nmp.sum(xFa[:,ja]) == 0.0:
            print('     Well! Seems that for variable '+L_VARO[jv]+', and algo '+calgo+', skin param has no impact')
            print('          ==> skipping difference plot...')

        else:

            fig = plt.figure(num = jv, figsize=size_fig, facecolor='w', edgecolor='k')
            ax1 = plt.axes([0.07, 0.22, 0.9, 0.75])

            ax1.set_xticks(vtime[::xticks_d])
            ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
            plt.xticks(rotation='60')

            plt.plot(vtime, xFa[:,ja], '-', color=l_color[ja], linestyle=l_style[ja], linewidth=l_width[ja], label=None, zorder=10+ja)

            ax1.set_ylim(-yrng,yrng) ; ax1.set_xlim(vtime[0],vtime[nbr-1])
            plt.ylabel(L_VARL[jv]+' ['+L_VUNT[jv]+']')

            ax1.grid(color='k', linestyle='-', linewidth=0.3)
            #plt.legend(bbox_to_anchor=(0.45, 0.2), ncol=1, shadow=True, fancybox=True)
            ax1.annotate(cvar_lnm+' ('+ctest+')', xy=(0.3, 0.97), xycoords='axes fraction',  bbox={'facecolor':'w', 'alpha':1., 'pad':10}, zorder=50, **font_inf)
            plt.savefig('diff_skin-noskin_'+L_VARO[jv]+'_'+calgo+'_'+ctest+'.'+fig_ext, dpi=int(rDPI), transparent=False)
            plt.close(jv)
