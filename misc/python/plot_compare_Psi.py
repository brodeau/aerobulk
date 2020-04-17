#!/usr/bin/env python
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-

# AeroBulk, 2020, L. Brodeau
#  (https://github.com/brodeau/aerobulk/)
#

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


dir_figs='.'
size_fig=(8,8)
fig_ext='png'

clr_red = '#AD0000'
clr_sat = '#ffed00'
clr_mod = '#008ab8'

rDPI=120.

L_ALGOS = [  'COARE'   , 'ECMWF'   , 'NCAR' , 'ANDREAS' ]
l_color = [  '#ffed00' , '#008ab8' , '0.6'  , '#AD0000' ] ; # colors to differentiate algos on the plot
l_width = [     5      ,    3      ,  2     ,  2        ] ; # line-width to differentiate algos on the plot
l_style = [    '-'     ,   '-'     , '--'   , '--'      ] ; # line-style
nb_algos = len(L_ALGOS) ; print(nb_algos)



# Getting arguments:
narg = len(sys.argv)
if narg != 2:
    print 'Usage: '+sys.argv[0]+' <output_from_test_psi_stab.x>.nc'; sys.exit(0)
cf_in = sys.argv[1]



if not path.exists(cf_in): print('ERROR: File '+cf_in+' does not exist !!!') ; sys.exit(0)

# Reading Psi profiles:
id_in = Dataset(cf_in)
vzeta = id_in.variables['zeta'][:]
#
nzeta = len(vzeta)
XPSI_M  = nmp.zeros( (nzeta,nb_algos) )
XPSI_H  = nmp.zeros( (nzeta,nb_algos) )
for ja in range(nb_algos):
    cv_psi_m = 'Psi_m_'+L_ALGOS[ja]
    cv_psi_h = 'Psi_h_'+L_ALGOS[ja]
    XPSI_M[:,ja] = id_in.variables[cv_psi_m][:]
    XPSI_H[:,ja] = id_in.variables[cv_psi_h][:]
id_in.close()



params = { 'font.family': 'Open Sans',
           'font.size':   40,
           'legend.fontsize': 13,
           'xtick.labelsize': 14,
           'ytick.labelsize': 14,
           'axes.labelsize':  16}
mpl.rcParams.update(params)
font_inf = { 'fontname':'Open Sans', 'fontweight':'normal', 'fontsize':18 }


fig = plt.figure(num=1, figsize=size_fig, facecolor='w', edgecolor='k')

ax1 = plt.axes([0.11, 0.09, 0.86, 0.88])

for ja in range(nb_algos):
    plt.plot(vzeta, XPSI_M[:,ja], '-', color=l_color[ja], linestyle=l_style[ja], linewidth=l_width[ja], label=L_ALGOS[ja], zorder=10+ja)

ax1.set_ylim(-30.,5.)
ax1.set_xlim(vzeta[0],vzeta[nzeta-1])
plt.ylabel(r'$\Psi_m(\zeta)$')
plt.xlabel(r'$\zeta$')

ax1.grid(color='k', linestyle='-', linewidth=0.3)
plt.legend(bbox_to_anchor=(0.45, 0.2), ncol=1, shadow=True, fancybox=True)
#ax1.annotate(cvar_lnm+' ('+ctest+')', xy=(0.3, 0.97), xycoords='axes fraction',  bbox={'facecolor':'w', 'alpha':1., 'pad':10}, zorder=50, **font_inf)
plt.savefig('Comparaison_Psi_m.'+fig_ext, dpi=int(rDPI), transparent=False)
plt.close(1)



sys.exit(0)


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
