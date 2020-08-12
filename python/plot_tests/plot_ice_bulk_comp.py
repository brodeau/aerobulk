#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-

# Post-diagnostic of test_aerobulk_buoy_series_ice.x

import sys
from os import path as path
import math
import numpy as nmp
from netCDF4 import Dataset,num2date
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

import clprn_tool as clt



dir_figs='.'
size_fig=(13,8)
size_fig0=(12,10)
fig_ext='png'

clr_red = '#AD0000'
clr_blu = '#3749A3'
clr_gre = '#548F64'
clr_sat = '#ffed00'
clr_mod = '#008ab8'

rDPI=100.

L_ALGOS = [ 'nemo'     ,     'an05'   ,        'lu12'        ,        'lg15'            ]
l_color = [   '0.1'    ,   clr_gre    ,       clr_blu        ,        clr_red           ] ; # colors to differentiate algos on the plot
l_width = [     1      ,       1      ,           1          ,          0.8             ] ; # line-width to differentiate algos on the plot
l_style = [    '-'     ,      '-'     ,          '-'         ,          '-'             ] ; # line-style
l_lgnm  = [ 'NEMO def.','Andreas 2005','Lupkes et al. (2012)','Lupkes & Gryanik (2015)' ]


L_VNEM  = [   'Qlat'    ,    'Qsen'     ,  'Tau'      ,   'Qlw'    ,   'z0'     ,   'Cd_i'    ,   'Rib_zu' ,   'CdN'     ,   'A'  ]
L_VARO  = [   'Qlat'    ,    'Qsen'     ,  'Tau'      ,   'Qlw'    ,   'z0'     ,   'Cd_i'    ,   'Rib'    ,   'CdN'     ,   'A'  ] ; # name of variable on figure
L_VARL  = [ r'$Q_{lat}$', r'$Q_{sens}$' , r'$|\tau|$' , r'$Q_{lw}$', r'$z_{0}$' , r'$C_{D}$'  , r'$Ri_{B}$', r'C_{D}^{N}',   'A'  ] ; # name of variable in latex mode
L_VUNT  = [ r'$W/m^2$'  , r'$W/m^2$'    , r'$N/m^2$'  , r'$W/m^2$' ,   r'$m$'   ,    r''      ,    r''     ,    r''      ,    ''  ]
L_VMAX  = [      50.     ,     50.       ,    0.5     ,     20.    ,     0.006   ,    3.5     ,    5.      ,    2.5      ,    1.  ]
L_VMIN  = [     -50.     ,    -50.       ,    0.      ,   -100.    ,     0.     ,     0.      ,   -5.      ,    1.       ,    0.  ]
L_ANOM  = [   True      ,    True       ,   True      ,    True    ,    True    ,    True     ,    False   ,    False    , False  ]


nb_algos = len(L_ALGOS)

# Getting arguments:
narg = len(sys.argv)
if narg != 2:
    print('Usage: '+sys.argv[0]+' <DIR_OUT_SASF>'); sys.exit(0)
cdir_data = sys.argv[1]



# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Populating and checking existence of files to be read
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def chck4f(cf):
    cmesg = 'ERROR: File '+cf+' does not exist !!!'
    if not path.exists(cf): print(cmesg) ; sys.exit(0)

cf_in = []
for ja in range(nb_algos):
    cfi = cdir_data+'/aerobulk_test_ice_series_'+L_ALGOS[ja]+'.nc'
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
vtime = vtime.astype(dtype='datetime64[D]')

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

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fig = plt.figure(num = jv, figsize=size_fig, facecolor='w', edgecolor='k')
    ax1 = plt.axes([0.08, 0.25, 0.9, 0.7])
    ax1.set_xticks(vtime[::xticks_d])
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
    plt.xticks(rotation='60', **font_x)

    for ja in range(nb_algos):
        plt.plot(vtime, xF[:,ja], '-', color=l_color[ja], linestyle=l_style[ja], \
                 linewidth=l_width[ja], label=l_lgnm[ja], zorder=10+ja)

    ax1.set_ylim(L_VMIN[jv], L_VMAX[jv]) ; ax1.set_xlim(vtime[0],vtime[Nt-1])
    plt.ylabel(L_VARL[jv]+' ['+L_VUNT[jv]+']')

    if L_VARO[jv] == 'Rib':
        plt.yscale('symlog', linthreshy=0.015)
    
    ax1.grid(color='k', linestyle='-', linewidth=0.3)
    plt.legend(loc='best', ncol=1, shadow=True, fancybox=True)
    ax1.annotate(cvar_lnm+' over sea-ice', xy=(0.3, 0.97), xycoords='axes fraction',  \
                 bbox={'facecolor':'w', 'alpha':1., 'pad':10}, zorder=50, **font_inf)
    plt.savefig(L_VARO[jv]+'.'+fig_ext, dpi=int(rDPI), transparent=False)
    plt.close(jv)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    if L_ANOM[jv]:

        for ja in range(nb_algos): xFa[:,ja] = xF[:,ja] - nmp.mean(xF,axis=1)

        if nmp.sum(nmp.abs(xFa[:,:])) == 0.0:
            print('     Well! Seems that for variable '+L_VARO[jv]+', choice of algo has no impact a all!')
            print('          ==> skipping anomaly plot...')

        else:

            yrng = clt.symetric_range( nmp.min(xFa) , nmp.max(xFa) )

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            fig = plt.figure(num = 10+jv, figsize=size_fig, facecolor='w', edgecolor='k')
            ax1 = plt.axes([0.08, 0.25, 0.9, 0.7])
            ax1.set_xticks(vtime[::xticks_d])
            ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
            plt.xticks(rotation='60', **font_x)

            for ja in range(nb_algos):
                plt.plot(vtime, xFa[:,ja], '-', color=l_color[ja], linewidth=l_width[ja], \
                         label=L_ALGOS[ja], zorder=10+ja)

            ax1.set_ylim(-yrng,yrng) ; ax1.set_xlim(vtime[0],vtime[Nt-1])
            plt.ylabel(L_VARL[jv]+' ['+L_VUNT[jv]+']')
            ax1.grid(color='k', linestyle='-', linewidth=0.3)
            plt.legend(bbox_to_anchor=(0.45, 0.2), ncol=1, shadow=True, fancybox=True)
            ax1.annotate('Anomaly of '+cvar_lnm, xy=(0.3, 0.97), xycoords='axes fraction',  \
                         bbox={'facecolor':'w', 'alpha':1., 'pad':10}, zorder=50, **font_inf)
            plt.savefig(L_VARO[jv]+'_anomaly.'+fig_ext, dpi=int(rDPI), transparent=False)
            plt.close(10+jv)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




for ja in range(nb_algos):

    i_cdn=7
    i_rib=6
    i_cd=5
    i_a=8

    # RiB
    id_in = Dataset(cf_in[ja])
    vrib = id_in.variables[L_VNEM[i_rib]][:] # only the center point of the 3x3 spatial domain!
    crib_lnm = id_in.variables[L_VNEM[i_rib]].long_name
    id_in.close()

    # Cd
    id_in = Dataset(cf_in[ja])
    vcd = id_in.variables[L_VNEM[i_cd]][:] # only the center point of the 3x3 spatial domain!
    ccd_lnm = id_in.variables[L_VNEM[i_cd]].long_name
    id_in.close()

    # Cdn
    id_in = Dataset(cf_in[ja])
    vcdn = id_in.variables[L_VNEM[i_cdn]][:] # only the center point of the 3x3 spatial domain!
    ccdn_lnm = id_in.variables[L_VNEM[i_cdn]].long_name
    id_in.close()

    # A (ice fraction)
    id_in = Dataset(cf_in[ja])
    va = id_in.variables[L_VNEM[i_a]][:] # only the center point of the 3x3 spatial domain!
    ca_lnm = id_in.variables[L_VNEM[i_a]].long_name
    id_in.close()


    rDPI=100.

    # CD(RiB)
    fig = plt.figure(num = jv, figsize=size_fig, facecolor='w', edgecolor='k')    
    ax1 = plt.axes([0.07, 0.1, 0.9, 0.83])
    plt.scatter(vrib, vcd, s=4, c='k', alpha=0.7)
    ax1.set_ylim(0.,3.) ; # Range for CD*1000
    plt.ylabel(r'$C_D$')
    ax1.set_xlim(-2.,2.) ; # Range for RiB
    plt.xlabel(r'$Ri_B$')
    ax1.grid(color='k', linestyle='-', linewidth=0.3)
    plt.savefig('scat_Cd_Rib_'+L_ALGOS[ja]+'.'+fig_ext, dpi=int(rDPI), transparent=False)
    plt.close(jv)


    # CD/CDN(RiB)
    fig = plt.figure(num = jv, figsize=size_fig, facecolor='w', edgecolor='k')    
    ax1 = plt.axes([0.07, 0.1, 0.9, 0.83])
    plt.scatter(vrib, vcd/vcdn, s=4, c='k', alpha=0.7)
    ax1.set_ylim(0.,2.) ; # Range for CD*1000
    plt.ylabel(r'$C_D/C_D^N$')
    ax1.set_xlim(-1.,1.) ; # Range for RiB
    plt.xlabel(r'$Ri_B$')
    ax1.grid(color='k', linestyle='-', linewidth=0.3)
    plt.savefig('scat_CdoCdN_Rib_'+L_ALGOS[ja]+'.'+fig_ext, dpi=int(rDPI), transparent=False)
    plt.close(jv)
    
    
    # CDN10 ( A )
    fig = plt.figure(num = jv, figsize=size_fig, facecolor='w', edgecolor='k')    
    ax1 = plt.axes([0.07, 0.1, 0.9, 0.83])
    plt.scatter(va, vcdn, s=4, c='k', alpha=0.7)
    ax1.set_ylim(1.,2.5) ; # Range for CDN*1000
    plt.ylabel(r'$C_D^{N10}$')
    ax1.set_xlim(0.,1.) ; # Range for A
    plt.xlabel(r'Ice concentration')
    ax1.grid(color='k', linestyle='-', linewidth=0.3)
    plt.savefig('scat_CdN_A_'+L_ALGOS[ja]+'.'+fig_ext, dpi=int(rDPI), transparent=False)
    plt.close(jv)
    
    
    
    
