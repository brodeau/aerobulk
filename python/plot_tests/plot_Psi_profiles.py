#!/usr/bin/env python
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-

# AeroBulk, 2020, L. Brodeau
#  (https://github.com/brodeau/aerobulk/)
#

import sys
from os import path as path
import math
import numpy as nmp
from netCDF4 import Dataset
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

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
font_axes = { 'fontname':'Arial',     'fontweight':'normal', 'fontsize':16 }


fig = plt.figure(num=1, figsize=size_fig, facecolor='w', edgecolor='k')
ax1 = plt.axes([0.11, 0.09, 0.86, 0.88])
#
for ja in range(nb_algos):
    plt.plot(vzeta, XPSI_M[:,ja], '-', color=l_color[ja], linestyle=l_style[ja], linewidth=l_width[ja], label=L_ALGOS[ja], zorder=10+ja)
ax1.set_ylim(-30.,5.)
ax1.set_xlim(vzeta[0],vzeta[nzeta-1])
plt.ylabel(r'$\Psi_m(\zeta)$', **font_axes)
plt.xlabel(r'$\zeta$', **font_axes)
ax1.grid(color='k', linestyle='-', linewidth=0.3)
plt.legend(bbox_to_anchor=(0.45, 0.2), ncol=1, shadow=True, fancybox=True)
#
plt.savefig('Comparaison_Psi_m.'+fig_ext, dpi=int(rDPI), transparent=False)
plt.close(1)



fig = plt.figure(num=1, figsize=size_fig, facecolor='w', edgecolor='k')
ax1 = plt.axes([0.11, 0.09, 0.86, 0.88])
#
for ja in range(nb_algos):
    plt.plot(vzeta, XPSI_H[:,ja], '-', color=l_color[ja], linestyle=l_style[ja], linewidth=l_width[ja], label=L_ALGOS[ja], zorder=10+ja)
ax1.set_ylim(-30.,5.)
ax1.set_xlim(vzeta[0],vzeta[nzeta-1])
plt.ylabel(r'$\Psi_h(\zeta)$', **font_axes)
plt.xlabel(r'$\zeta$', **font_axes)
ax1.grid(color='k', linestyle='-', linewidth=0.3)
plt.legend(bbox_to_anchor=(0.45, 0.2), ncol=1, shadow=True, fancybox=True)
#
plt.savefig('Comparaison_Psi_h.'+fig_ext, dpi=int(rDPI), transparent=False)
plt.close(1)
