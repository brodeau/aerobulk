#!/usr/bin/env python
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-

# Post-diagnostic of test_aerobulk_buoy_series_ice.x

import sys
from os import path as path
#from string import replace
import math
import numpy as nmp
from netCDF4 import Dataset,num2date
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

import clprn_tool as clt

reload(sys)
sys.setdefaultencoding('utf8')



dir_figs='.'
size_fig=(13,7)
fig_ext='png'

clr_red = '#AD0000'
clr_blu = '#3749A3'
clr_gre = '#548F64'
clr_sat = '#ffed00'
clr_mod = '#008ab8'

#rDPI=400.
rDPI=150.

L_ALGOS = [ 'nemo'     ,     'an05'     ,    'lu15'     ,  'best'  ]
l_color = [   '0.1'    ,   clr_gre      , clr_blu       ,  clr_red ] ; # colors to differentiate algos on the plot
l_width = [     1      ,       1        ,     0.6       ,    0.5   ] ; # line-width to differentiate algos on the plot
l_style = [    '-'     ,      '-'       ,     '-'       ,  '--'    ] ; # line-style
l_lgnm  = [ 'NEMO def.','Andreas (2005)','Lupkes (2015)','Brodeau' ]

#L_VNEM  = [   'qla'     ,     'qsb'     ,     'qt'     ,   'qlw'     ,  'taum'     ,    'dt_skin'         ]
#L_VARO  = [   'Qlat'    ,    'Qsen'     ,     'Qnet'   ,   'Qlw'     ,  'Tau'      ,    'dT_skin'         ] ; # name of variable on figure
#L_VARL  = [ r'$Q_{lat}$', r'$Q_{sens}$' , r'$Q_{net}$' , r'$Q_{lw}$' , r'$|\tau|$' , r'$\Delta T_{skin}$' ] ; # name of variable in latex mode
#L_VUNT  = [ r'$W/m^2$'  , r'$W/m^2$'    , r'$W/m^2$'   , r'$W/m^2$'  , r'$N/m^2$'  ,      'K'             ]
#L_VMAX  = [     75.     ,     25.       ,    800.      ,     25.     ,    1.2      ,      -0.7            ]
#L_VMIN  = [   -250.     ,    -25.       ,   -400.      ,   -150.     ,    0.       ,       0.7            ]
#L_ANOM  = [   True      ,    True       ,    True      ,    True     ,   True      ,      False           ]

L_VNEM  = [   'Qlat'    ,    'Qsen'     ,  'Tau'      ,   'Qlw'    ,   'z0'     ,   'Cd'    ,   'Rib'    ]
L_VARO  = [   'Qlat'    ,    'Qsen'     ,  'Tau'      ,   'Qlw'    ,   'z0'     ,   'Cd'    ,   'Rib'    ] ; # name of variable on figure
L_VARL  = [ r'$Q_{lat}$', r'$Q_{sens}$' , r'$|\tau|$' , r'$Q_{lw}$', r'$z_{0}$' , r'$C_{D}$', r'$Ri_{B}$'] ; # name of variable in latex mode
L_VUNT  = [ r'$W/m^2$'  , r'$W/m^2$'    , r'$N/m^2$'  , r'$W/m^2$' ,   r'$m$'   ,    r''    ,    r''     ]
L_VMAX  = [      50.     ,     50.       ,    0.5     ,     20.    ,     0.01   ,    3.5    ,    5.      ]
L_VMIN  = [     -50.     ,    -50.       ,    0.      ,   -100.    ,     0.     ,     0.    ,   -5.      ]
L_ANOM  = [   True      ,    True       ,   True      ,    True    ,    True    ,    True   ,    False   ]


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

cf_in = []
for ja in range(nb_algos):
    cfi = cdir_data+'/lolo_'+L_ALGOS[ja]+'.nc'
    chck4f(cfi)
    cf_in.append(cfi)
print('Files we are goin to use:')
for ja in range(nb_algos): print(cf_in[ja])
#-----------------------------------------------------------------



# Getting time array from the first file:
id_in = Dataset(cf_in[0])
vt = id_in.variables['time'][:]
cunit_t = id_in.variables['time'].units
clndr_t = id_in.variables['time'].calendar
id_in.close()
Nt = len(vt)
print(' "time" => units = '+cunit_t+', calendar = "'+clndr_t+'"')

vtime = num2date(vt, units=cunit_t) ; # something understandable!

ii=Nt/300
ib=max(ii-ii%10,1)
xticks_d=int(30*ib)

font_inf = { 'fontname':'Open Sans', 'fontweight':'normal', 'fontsize':14 }

nb_var = len(L_VNEM)

xF  = nmp.zeros((Nt,nb_algos))
xFa = nmp.zeros((Nt,nb_algos))




for jv in range(nb_var):
    print('\n *** Treating variable: '+L_VARO[jv]+' !')

    for ja in range(nb_algos):
        #
        id_in = Dataset(cf_in[ja])
        xF[:,ja] = id_in.variables[L_VNEM[jv]][:] # only the center point of the 3x3 spatial domain!
        if ja == 0: cvar_lnm = id_in.variables[L_VNEM[jv]].long_name
        id_in.close()

    fig = plt.figure(num = jv, figsize=size_fig, facecolor='w', edgecolor='k')

    ax1 = plt.axes([0.07, 0.22, 0.9, 0.75])

    ax1.set_xticks(vtime[::xticks_d])
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
    plt.xticks(rotation='60')

    for ja in range(nb_algos):
        plt.plot(vtime, xF[:,ja], '-', color=l_color[ja], linestyle=l_style[ja], linewidth=l_width[ja], label=l_lgnm[ja], zorder=10+ja)

    ax1.set_ylim(L_VMIN[jv], L_VMAX[jv]) ; ax1.set_xlim(vtime[0],vtime[Nt-1])
    plt.ylabel(L_VARL[jv]+' ['+L_VUNT[jv]+']')

    ax1.grid(color='k', linestyle='-', linewidth=0.3)
    plt.legend(bbox_to_anchor=(0.45, 0.2), ncol=1, shadow=True, fancybox=True)
    ax1.annotate(cvar_lnm+' over sea-ice', xy=(0.3, 0.97), xycoords='axes fraction',  bbox={'facecolor':'w', 'alpha':1., 'pad':10}, zorder=50, **font_inf)
    plt.savefig(L_VARO[jv]+'.'+fig_ext, dpi=int(rDPI), transparent=False)
    plt.close(jv)



    if L_ANOM[jv]:

        for ja in range(nb_algos): xFa[:,ja] = xF[:,ja] - nmp.mean(xF,axis=1)

        if nmp.sum(nmp.abs(xFa[:,:])) == 0.0:
            print('     Well! Seems that for variable '+L_VARO[jv]+', choice of algo has no impact a all!')
            print('          ==> skipping anomaly plot...')
        
        else:

            yrng = clt.symetric_range( nmp.min(xFa) , nmp.max(xFa) )

            fig = plt.figure(num = 10+jv, figsize=size_fig, facecolor='w', edgecolor='k')
            ax1 = plt.axes([0.07, 0.22, 0.9, 0.75])

            ax1.set_xticks(vtime[::xticks_d])
            ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
            plt.xticks(rotation='60')

            for ja in range(nb_algos):
                plt.plot(vtime, xFa[:,ja], '-', color=l_color[ja], linewidth=l_width[ja], label=L_ALGOS[ja], zorder=10+ja)

            ax1.set_ylim(-yrng,yrng) ; ax1.set_xlim(vtime[0],vtime[Nt-1])
            plt.ylabel(L_VARL[jv]+' ['+L_VUNT[jv]+']')
            ax1.grid(color='k', linestyle='-', linewidth=0.3)
            plt.legend(bbox_to_anchor=(0.45, 0.2), ncol=1, shadow=True, fancybox=True)
            ax1.annotate('Anomaly of '+cvar_lnm, xy=(0.3, 0.97), xycoords='axes fraction',  bbox={'facecolor':'w', 'alpha':1., 'pad':10}, zorder=50, **font_inf)
            plt.savefig(L_VARO[jv]+'_anomaly.'+fig_ext, dpi=int(rDPI), transparent=False)
            plt.close(10+jv)









for ja in range(nb_algos):

    i_rib=6
    i_cd=5

    # Absciss = RiB
    id_in = Dataset(cf_in[ja])
    vx = id_in.variables[L_VNEM[i_rib]][:] # only the center point of the 3x3 spatial domain!
    crib_lnm = id_in.variables[L_VNEM[i_rib]].long_name
    id_in.close()

    # Ordonates = Cd
    id_in = Dataset(cf_in[ja])
    vcd = id_in.variables[L_VNEM[i_cd]][:] # only the center point of the 3x3 spatial domain!
    ccd_lnm = id_in.variables[L_VNEM[i_cd]].long_name
    id_in.close()



    fig = plt.figure(num = jv, figsize=size_fig, facecolor='w', edgecolor='k')    
    ax1 = plt.axes([0.07, 0.22, 0.9, 0.75])
    
    #ax1.set_xticks(vtime[::xticks_d])
    #ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
    #plt.xticks(rotation='60')
    
    plt.scatter(vx, vcd, s=4, c='k', alpha=0.7)
    #label=l_lgnm[ja], zorder=10+ja)
    
    ax1.set_ylim(0.,6.) ; # Range for CD*1000
    plt.ylabel(r'$C_D$')
    ax1.set_xlim(-2.,2.) ; # Range for RiB
    plt.xlabel(r'$Ri_B$')
    
    ax1.grid(color='k', linestyle='-', linewidth=0.3)
    #plt.legend(bbox_to_anchor=(0.45, 0.2), ncol=1, shadow=True, fancybox=True)
    #ax1.annotate(cvar_lnm+' over sea-ice', xy=(0.3, 0.97), xycoords='axes fraction',  \
    #    bbox={'facecolor':'w', 'alpha':1., 'pad':10}, zorder=50, **font_inf)
    plt.savefig('Cd_Rib_'+L_ALGOS[ja]+'.'+fig_ext, dpi=int(rDPI), transparent=False)
    plt.close(jv)
    
    
    
    
