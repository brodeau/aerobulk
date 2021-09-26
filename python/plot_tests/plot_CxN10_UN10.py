#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
# AeroBulk, 2020, L. Brodeau
#  (https://github.com/brodeau/aerobulk/)
#
# Post-diagnostic of test_cx_vs_wind.x

import sys
from os import getenv,mkdir,path
import numpy as nmp
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

l_multi_fig = False

#fig_ext = 'png' ; lshow_xylabs = True
fig_ext = 'svg' ; lshow_xylabs = True

AEROBULK_HOME = getenv('AEROBULK_HOME')
if AEROBULK_HOME is None:
    print('The environment variable "AEROBULK_HOME" must be set!')
    sys.exit(0)
cdir_in = AEROBULK_HOME+'/dat'    
print('\n *** Will look for data into '+AEROBULK_HOME+' !')

# Directory where to save figures:
fig_dir = AEROBULK_HOME+'/figures'
if not path.isdir(fig_dir):
    try:
        mkdir(fig_dir)
    except OSError:
        print ("Creation of the directory %s failed" % fig_dir)
    else:
        print ("Successfully created the directory %s " % fig_dir)


#vtv_u = [ '-0050','-0300','-1000' ]
#vtv_s = [ '+0050','+0300','+1000' ]



#valgo_nm    = [ 'ncar'   , 'coare3p0' , 'coare3p6' , 'ecmwf'  , 'andreas' ]
#valgo_DN    = [ 'NCAR'   , 'COARE 3.0', 'COARE 3.6', 'ECMWF'  , 'Andreas' ]
#vcolor      = [ '#3465a4', '#cc0000'  ,     'k'    , '#58FAAC',   'orange'  ]

#valgo_nm    = [ 'ncar'   , 'coare3p0' , 'coare3p6', 'ecmwf'   , 'andreas' ]
#valgo_DN    = [ 'NCAR'   , 'COARE 3.0', 'COARE'   ,  'ECMWF'  , 'Andreas' ]
#vcolor      = [ '#3465a4', '#D9344F'  ,    'k'    , '#58FAAC' , '#cc0000' ]
#vlines      = [   '-'    ,    '--'    ,    '-'    ,    '-'    ,   '--'    ]

valgo_nm    = [ 'ncar'   ,  'ecmwf' , 'coare3p0' , 'coare3p6'  ]
valgo_DN    = [ 'NCAR'   ,  'ECMWF' , 'COARE 3.0', 'COARE 3.6' ]
#vcolor     = [ '#3465a4', '#58FAAC', '#D9344F'  ,    'k'      ]
vcolor =      [  '0.6'   , '#008ab8', '#ffed00'  , '#C46637'   ]
vlines      = [   '-'    ,    '-'   ,    '--'    ,    '-'      ]

#valgo_nm    = [ 'ncar'   , 'coare3p6' ,  'ecmwf'  ]
#valgo_DN    = [ 'NCAR'   , 'COARE 3.6',  'ECMWF'  ]
#vcolor      = [ '#3465a4',    'k'     , '#58FAAC' ]
#vlines      = [   '-'    ,    '-'     ,    '-'    ]


nb_algo = len(valgo_nm)



############################################################################
def read_ascii_column(cfile, ivcol2read):
    if not path.exists(cfile):
        print('\n  ERROR: file '+cfile+' is missing !')
        sys.exit(0)
    f = open(cfile, 'r')
    cread_lines = f.readlines()
    f.close()
    nbcol = len(ivcol2read)
    jl = 0
    for ll in cread_lines:
        ls = ll.split()
        if ls[0] != '#': jl = jl + 1
    nbl = jl
    Xout  = nmp.zeros((nbcol,nbl))
    jl = -1
    for ll in cread_lines:
        ls = ll.split()
        if ls[0] != '#':
            jl = jl+1
            jc = -1
            for icol in ivcol2read:
                jc = jc+1
                Xout[jc,jl] = float(ls[icol])
    return Xout
############################################################################


cf_in = cdir_in+'/Neutral_coeff_U10N_'+valgo_nm[0]+'.dat'
xdum = read_ascii_column(cf_in, [0,1])
print('\n\n\n')


nU = len(xdum[0,:])

print(' *** nU =', nU)


xcd_u = nmp.zeros((2,nU,nb_algo))
for ja in range(nb_algo):
    xcd_u[:,:,ja] = read_ascii_column(cdir_in+'/Neutral_coeff_U10N_'+valgo_nm[ja]+'.dat', [0,1])

xch_u = nmp.zeros((2,nU,nb_algo))
for ja in range(nb_algo):
    xch_u[:,:,ja] = read_ascii_column(cdir_in+'/Neutral_coeff_U10N_'+valgo_nm[ja]+'.dat', [0,2])

xce_u = nmp.zeros((2,nU,nb_algo))
for ja in range(nb_algo):
    xce_u[:,:,ja] = read_ascii_column(cdir_in+'/Neutral_coeff_U10N_'+valgo_nm[ja]+'.dat', [0,3])


# x/y-axis range and increment:

vu_rng   = [ 0., 35., 1. ]
cd_u_rng = [0.8, 3. , 0.2]
ch_u_rng = [0.8, 1.8, 0.2]
ce_u_rng = [0.8, 1.8, 0.2]



DPI_FIG = 100

rscale  = 1.25
params = { 'font.family': 'Open Sans',
           'font.size':   40,
           'legend.fontsize': rscale*15,
           'xtick.labelsize': rscale*14,
           'ytick.labelsize': rscale*14,
           'axes.labelsize':  rscale*16}
mpl.rcParams.update(params)
#
font_ttl = { 'fontname':'Open Sans', 'fontweight':'normal', 'fontsize':rscale*22 }
font_ylb = { 'fontname':'Open Sans', 'fontweight':'normal', 'fontsize':rscale*16 }
font_xlb = { 'fontname':'Open Sans', 'fontweight':'normal', 'fontsize':rscale*18 }
font_lab = { 'fontname':'Open Sans', 'fontweight':'normal', 'fontsize':rscale*22 }
#font_info= { 'fontname':'Open Sans', 'fontweight':'normal', 'fontsize':rscale*16 }
#font_lab = { 'fontname':'Open Sans', 'fontweight':'normal', 'fontsize':rscale*16 }


y_tit = 0.92
y_lbl = 0.9

rwidt = 0.425
rheig = 0.25
rleft = 0.063
rbott = 0.04
rbext = 0.09



vU = nmp.zeros(nU)
vU[:] = xcd_u[0,:,0]



print('vU =', vU[:])


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def plot_var( nba, cvar_nm, vun10, xCd, xCe, valg_dn, x_rng=None, y_rng=None, istab=0, cunit='' ):
    #
    fig = plt.figure(num=1, figsize=(14.,10.), facecolor='w', edgecolor='k')
    ax1 = plt.axes([0.08, 0.08, 0.9, 0.9])    
    ax1.grid(color='k', linestyle='-', linewidth=0.2, zorder=0.1)
    #ax1.annotate(cstblt, xy=(0.4, 1.06), xycoords='axes fraction', **font_ttl)
    for ja in range(nba):
        ax1.plot( vun10, xCd[1,:,ja], label=valg_dn[ja], \
                  linewidth=5., color=vcolor[ja], linestyle=vlines[ja], zorder=2. )
    for ja in range(nba):
        ax1.plot( vun10, xCe[1,:,ja], label=None, \
                  linewidth=2.5, color=vcolor[ja], linestyle=vlines[ja], zorder=0.75 )
        #
    if nmp.shape(x_rng) == (3,):
        vx_ticks = nmp.arange(x_rng[0], x_rng[1]+x_rng[2], x_rng[2])
        ax1.set_xlim(x_rng[0], x_rng[1])
    if nmp.shape(y_rng) == (3,):
        vy_ticks = nmp.arange(y_rng[0], y_rng[1]+y_rng[2], y_rng[2]) ; plt.yticks( vy_ticks )
        ax1.set_ylim(y_rng[0], y_rng[1])
    #
    ctit = r'$C_{D}^{N10}$, $C_{E}^{N10}$'
    plt.text( 0.5, y_tit, ctit, horizontalalignment='center', transform = ax1.transAxes, \
              bbox=dict(boxstyle="square", fc='w'), **font_ttl)
    if not cunit=='': plt.ylabel(r''+cunit, **font_ylb)
    plt.xlabel(r'$U_{N10}$ [m/s]',     **font_xlb)
    #plt.legend(bbox_to_anchor=(0.34, 0.16), ncol=3, shadow=False, fancybox=True)
    plt.legend(bbox_to_anchor=(0.15, 0.63), ncol=1, shadow=False, fancybox=True)
    #
    plt.savefig(fig_dir+'/'+cvar_nm+'.'+fig_ext, dpi=DPI_FIG, facecolor='w', edgecolor='w', orientation='portrait')
    plt.close(1)
    return 0
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# Wind speed is x-axis in following plots:

iplt = plot_var( nb_algo, 'CXN10_vs_UN10', vU, xcd_u, xce_u, valgo_DN, x_rng=vu_rng, y_rng=cd_u_rng, istab=0, cunit='$10^{-3}$' )

print('\n *** Check figure into '+fig_dir+'/ !!!\n')

