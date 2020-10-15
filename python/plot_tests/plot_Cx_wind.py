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
sst_ref = 15.

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

fig_ext = 'png' ; lshow_xylabs = True
#fig_ext = 'svg' ; lshow_xylabs = False

vtv_u = [ '-0050','-0300','-1000' ]
vtv_s = [ '+0050','+0300','+1000' ]
ntv = len(vtv_u)



#valgo_nm    = [ 'ncar'   , 'coare3p0' , 'coare3p6' , 'ecmwf'  , 'andreas' ]
#valgo_DN    = [ 'NCAR'   , 'COARE 3.0', 'COARE 3.6', 'ECMWF'  , 'Andreas' ]
#vcolor      = [ '#3465a4', '#cc0000'  ,     'k'    , '#58FAAC',   'orange'  ]

valgo_nm    = [ 'ncar'   , 'coare3p6', 'ecmwf'   , 'andreas' ]
valgo_DN    = [ 'NCAR'   , 'COARE'   ,  'ECMWF'  , 'Andreas' ]
vcolor      = [ '#3465a4',    'k'    , '#58FAAC' , '#cc0000' ]
vlines      = [   '-'    ,    '-'    ,    '-'    ,   '--'    ]

#valgo_nm    = [ 'ncar'   , 'coare3p6' ]
#valgo_DN    = [ 'NCAR'   , 'COARE'    ]
#vcolor      = [ '#3465a4',    'k'     ]


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


vrt_u = nmp.zeros(ntv) ; vrt_s = nmp.zeros(ntv)


csst = str(int(sst_ref))

cf_in = cdir_in+'/cd_dtv_-0500_sst_'+csst+'_'+valgo_nm[0]+'.dat'
xdum = read_ascii_column(cf_in, [0,1])
print('\n\n\n')


nU = len(xdum[0,:])

xcd_u = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    vrt_u[jtv] = float(vtv_u[jtv])/100. # dTv as a float
    for ja in range(nb_algo):
        xcd_u[:,:,jtv,ja] = read_ascii_column(cdir_in+'/cd_dtv_'+vtv_u[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])

xcd_s = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    vrt_s[jtv] = float(vtv_s[jtv])/100. # dTv as a float
    for ja in range(nb_algo):
        xcd_s[:,:,jtv,ja] = read_ascii_column(cdir_in+'/cd_dtv_'+vtv_s[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])

xch_u = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    for ja in range(nb_algo):
        xch_u[:,:,jtv,ja] = read_ascii_column(cdir_in+'/ch_dtv_'+vtv_u[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])

xch_s = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    for ja in range(nb_algo):
        xch_s[:,:,jtv,ja] = read_ascii_column(cdir_in+'/ch_dtv_'+vtv_s[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])


xce_u = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    for ja in range(nb_algo):
        xce_u[:,:,jtv,ja] = read_ascii_column(cdir_in+'/ce_dtv_'+vtv_u[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])

xce_s = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    for ja in range(nb_algo):
        xce_s[:,:,jtv,ja] = read_ascii_column(cdir_in+'/ce_dtv_'+vtv_s[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])


##### L ####
xlo_u = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    #vrt_u[jtv] = float(vtv_u[jtv])/100. # dTv as a float
    for ja in range(nb_algo):
        xlo_u[:,:,jtv,ja] = read_ascii_column(cdir_in+'/lo_dtv_'+vtv_u[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])

xlo_s = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    #vrt_s[jtv] = float(vtv_s[jtv])/100. # dTv as a float
    for ja in range(nb_algo):
        xlo_s[:,:,jtv,ja] = read_ascii_column(cdir_in+'/lo_dtv_'+vtv_s[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])

##### UN10 ####
xun_u = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    #vrt_u[jtv] = float(vtv_u[jtv])/100. # dTv as a float
    for ja in range(nb_algo):
        xun_u[:,:,jtv,ja] = read_ascii_column(cdir_in+'/un_dtv_'+vtv_u[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])
        xun_u[1,:,jtv,ja] = xun_u[1,:,jtv,ja] ; # => UN10

xun_s = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    #vrt_s[jtv] = float(vtv_s[jtv])/100. # dTv as a float
    for ja in range(nb_algo):
        xun_s[:,:,jtv,ja] = read_ascii_column(cdir_in+'/un_dtv_'+vtv_s[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])
        xun_s[1,:,jtv,ja] = xun_s[1,:,jtv,ja]; # => UN10

##### u* ####
xus_u = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    #vrt_u[jtv] = float(vtv_u[jtv])/100. # dTv as a float
    for ja in range(nb_algo):
        xus_u[:,:,jtv,ja] = read_ascii_column(cdir_in+'/us_dtv_'+vtv_u[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])

xus_s = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    #vrt_s[jtv] = float(vtv_s[jtv])/100. # dTv as a float
    for ja in range(nb_algo):
        xus_s[:,:,jtv,ja] = read_ascii_column(cdir_in+'/us_dtv_'+vtv_s[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])

##### Rib ####
xri_u = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    #vrt_u[jtv] = float(vtv_u[jtv])/100. # dTv as a float
    for ja in range(nb_algo):
        xri_u[:,:,jtv,ja] = read_ascii_column(cdir_in+'/ri_dtv_'+vtv_u[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])

xri_s = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    #vrt_s[jtv] = float(vtv_s[jtv])/100. # dTv as a float
    for ja in range(nb_algo):
        xri_s[:,:,jtv,ja] = read_ascii_column(cdir_in+'/ri_dtv_'+vtv_s[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])

##### z0 ####
xz0_u = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    #vrt_u[jtv] = float(vtv_u[jtv])/100. # dTv as a float
    for ja in range(nb_algo):
        xz0_u[:,:,jtv,ja] = read_ascii_column(cdir_in+'/z0_dtv_'+vtv_u[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])

xz0_s = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    #vrt_s[jtv] = float(vtv_s[jtv])/100. # dTv as a float
    for ja in range(nb_algo):
        xz0_s[:,:,jtv,ja] = read_ascii_column(cdir_in+'/z0_dtv_'+vtv_s[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])



# x/y-axis range and increment:

vu_rng   = [ 0., 34., 1. ]

cd_u_rng = [0.8, 3. , 0.2]
cd_s_rng = [0.0, 2.5, 0.5]

ch_u_rng = [0.8, 1.8, 0.2]
ch_s_rng = [0.0, 1.4, 0.2]

ce_u_rng = [0.8, 1.8, 0.2]
ce_s_rng = [0.0, 1.6, 0.2]


DPI_FIG = 100

rscale  = 1.25
params = { 'font.family': 'Open Sans',
           'font.size':   40,
           'legend.fontsize': rscale*10,
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
vU[:] = xcd_u[0,:,jtv,0]



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def plot_var_u( nvt, nba, cvar_nm, vu10, xF, vstab, valg_dn, x_rng=None, y_rng=None, istab=0, cunit='' ):
    #
    cstblt = 'unstable'
    csgn   = '-'
    if istab==1:
        cstblt = 'stable'
        csgn   = '+'
    #
    fig = plt.figure(num=1, figsize=(14.,10.), facecolor='w', edgecolor='k')
    ax1 = plt.axes([0.08, 0.08, 0.9, 0.9])    
    ax1.grid(color='k', linestyle='-', linewidth=0.2, zorder=0.1)
    ax1.annotate(cstblt, xy=(0.4, 1.06), xycoords='axes fraction', **font_ttl)
    for jtv in range(nvt):
        for ja in range(nba):
            ax1.plot( vu10, xF[1,:,jtv,ja], label=valg_dn[ja]+r' $\Delta\Theta=$'+csgn+str(vstab[jtv])+'K', \
                      linewidth=3.5-jtv*0.7, color=vcolor[ja], linestyle=vlines[ja], zorder=0.75 )
    #
    if nmp.shape(x_rng) == (3,):
        vx_ticks = nmp.arange(x_rng[0], x_rng[1]+x_rng[2], x_rng[2])
        ax1.set_xlim(x_rng[0], x_rng[1])
    if nmp.shape(y_rng) == (3,):
        vy_ticks = nmp.arange(y_rng[0], y_rng[1]+y_rng[2], y_rng[2]) ; plt.yticks( vy_ticks )
        ax1.set_ylim(y_rng[0], y_rng[1])
    #
    ctit = r'$'+cvar_nm+'$ --'+cstblt+'--'
    plt.text(0.5, y_tit, ctit, horizontalalignment='center', transform = ax1.transAxes, bbox=dict(boxstyle="square", fc='w'), **font_ttl)
    if not cunit=='': plt.ylabel(r''+cunit, **font_ylb)
    plt.xlabel(r'$U_{10m}$ [m/s]',     **font_xlb)
    plt.legend(bbox_to_anchor=(0.34, 0.16), ncol=3, shadow=False, fancybox=True)
    #
    plt.savefig(fig_dir+'/'+cvar_nm+'_'+cstblt+'.'+fig_ext, dpi=DPI_FIG, facecolor='w', edgecolor='w', orientation='portrait')
    plt.close(1)
    return 0
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def plot_scat( nba, cvar_nm, xF1, xF2, valg_dn, x_rng=None, y_rng=None, istab=0, cxunit='', cyunit='' ):
    #
    cstblt = 'unstable'
    csgn   = '-'
    if istab==1:
        cstblt = 'stable'
        csgn   = '+'
    #
    fig = plt.figure(num=1, figsize=(14.,10.), facecolor='w', edgecolor='k')
    ax1 = plt.axes([0.07, 0.08, 0.9, 0.9])    
    ax1.grid(color='k', linestyle='-', linewidth=0.2, zorder=0.1)
    ax1.annotate(cstblt, xy=(0.4, 1.06), xycoords='axes fraction', **font_ttl)
    #
    for ja in range(nba):
        plt.scatter( xF1[1,:,:,ja], xF2[1,:,:,ja], marker='.', label=valg_dn[ja] )
    
    if nmp.shape(x_rng) == (3,):
        vx_ticks = nmp.arange(x_rng[0], x_rng[1]+x_rng[2], x_rng[2]) ; plt.xticks( vx_ticks )
        ax1.set_xlim(x_rng[0], x_rng[1])
    #
    if nmp.shape(y_rng) == (3,):
        vy_ticks = nmp.arange(y_rng[0], y_rng[1]+y_rng[2], y_rng[2]) ; plt.yticks( vy_ticks )
        ax1.set_ylim(y_rng[0], y_rng[1])
    #
    ctit = r'$'+cvar_nm+'$ --'+cstblt+'--'
    plt.text(0.5, y_tit, ctit, horizontalalignment='center', transform = ax1.transAxes, bbox=dict(boxstyle="square", fc='w'), **font_ttl)
    plt.ylabel(r''+cyunit, **font_ylb)
    plt.xlabel(r''+cxunit, **font_xlb)
    plt.legend(bbox_to_anchor=(0.7, 0.16), ncol=1, shadow=False, fancybox=True)
    #
    plt.savefig(fig_dir+'/'+cvar_nm+'_'+cstblt+'.'+fig_ext, dpi=DPI_FIG, facecolor='w', edgecolor='w', orientation='portrait')
    plt.close(1)
    return 0
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







# Wind speed is x-axis in following plots:

iplt = plot_var_u( ntv, nb_algo, 'C_D', vU, xcd_u, vrt_u, valgo_DN, x_rng=vu_rng, y_rng=cd_u_rng, istab=0, cunit='$10^{-3}$' )
iplt = plot_var_u( ntv, nb_algo, 'C_D', vU, xcd_s, vrt_s, valgo_DN, x_rng=vu_rng, y_rng=cd_s_rng, istab=1, cunit='$10^{-3}$' )

iplt = plot_var_u( ntv, nb_algo, 'C_H', vU, xch_u, vrt_u, valgo_DN, x_rng=vu_rng, y_rng=ch_u_rng, istab=0, cunit='$10^{-3}$' )
iplt = plot_var_u( ntv, nb_algo, 'C_H', vU, xch_s, vrt_s, valgo_DN, x_rng=vu_rng, y_rng=ch_s_rng, istab=1, cunit='$10^{-3}$' )

iplt = plot_var_u( ntv, nb_algo, 'C_E', vU, xce_u, vrt_u, valgo_DN, x_rng=vu_rng, y_rng=ce_u_rng, istab=0, cunit='$10^{-3}$' )
iplt = plot_var_u( ntv, nb_algo, 'C_E', vU, xce_s, vrt_s, valgo_DN, x_rng=vu_rng, y_rng=ce_s_rng, istab=1, cunit='$10^{-3}$' )

iplt = plot_var_u( ntv, nb_algo, 'L',   vU, xlo_u, vrt_s, valgo_DN, x_rng=vu_rng, y_rng=[-500.,   0., 50.], istab=0, cunit='[m]' )
iplt = plot_var_u( ntv, nb_algo, 'L',   vU, xlo_s, vrt_s, valgo_DN, x_rng=vu_rng, y_rng=[   0., 500., 50.], istab=1, cunit='[m]' )

iplt = plot_var_u( ntv, nb_algo, 'u*',  vU, xus_u, vrt_s, valgo_DN, x_rng=vu_rng, y_rng=[0., 2., 0.1], istab=0, cunit='[m/s]' )
iplt = plot_var_u( ntv, nb_algo, 'u*',  vU, xus_s, vrt_s, valgo_DN, x_rng=vu_rng, y_rng=[0., 2., 0.1], istab=1, cunit='[m/s]' )

iplt = plot_var_u( ntv, nb_algo, 'Ri_B',   vU, xri_u, vrt_s, valgo_DN, x_rng=vu_rng, y_rng=[-10.,  0., 1], istab=0, cunit='' )
iplt = plot_var_u( ntv, nb_algo, 'Ri_B',   vU, xri_s, vrt_s, valgo_DN, x_rng=vu_rng, y_rng=[  0., 10., 1], istab=1, cunit='' )

iplt = plot_var_u( ntv, nb_algo, 'z_0',   vU, xz0_u, vrt_s, valgo_DN, x_rng=vu_rng, y_rng=[0., 0.005, 0.001], istab=0, cunit='[m]' )
iplt = plot_var_u( ntv, nb_algo, 'z_0',   vU, xz0_s, vrt_s, valgo_DN, x_rng=vu_rng, y_rng=[0., 0.005, 0.001], istab=1, cunit='[m]' )





### x-axis is now neutral-stability wind speed:

# u*(UN10)    (must be more or less identical for stable and unstable !!!, but good to keep boths for debugging purposes...)
iplt = plot_scat( nb_algo, 'u*(UN10)', xun_u, xus_u, valgo_DN, x_rng=[0,30,1], y_rng=[0,1.5,0.1], istab=0, cxunit='$U_{N10}$ [m/s]', cyunit='u* [m/s]' )
iplt = plot_scat( nb_algo, 'u*(UN10)', xun_s, xus_s, valgo_DN, x_rng=[0,30,1], y_rng=[0,1.5,0.1], istab=1, cxunit='$U_{N10}$ [m/s]', cyunit='u* [m/s]' )















print('\n *** Check figures into '+fig_dir+'/ !!!\n')

