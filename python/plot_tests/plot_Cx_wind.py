#!/usr/bin/env python
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-

# Post-diagnostic of test_cx_vs_wind.x

import sys
from os import getenv
import numpy as nmp
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import clprn_tool as clt

l_multi_fig = False

sst_ref = 15.


AEROBULK_HOME = getenv('AEROBULK_HOME')
if AEROBULK_HOME is None:
    print('The environment variable "AEROBULK_HOME" must be set!')
    sys.exit(0)
cdir_in = AEROBULK_HOME+'/dat'    
print('\n *** Will look for data into '+AEROBULK_HOME+' !')



fig_frmt = 'png' ; lshow_xylabs = True
#fig_frmt = 'svg' ; lshow_xylabs = False

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


vrt_u = nmp.zeros(ntv) ; vrt_s = nmp.zeros(ntv)


csst = str(int(sst_ref))

cf_in = cdir_in+'/cd_dtv_-0500_sst_'+csst+'_'+valgo_nm[0]+'.dat'
xdum = clt.read_ascii_column(cf_in, [0,1])
print '\n\n\n'


nU = len(xdum[0,:])

xcd_u = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    vrt_u[jtv] = float(vtv_u[jtv])/100. # dTv as a float
    for ja in range(nb_algo):
        xcd_u[:,:,jtv,ja] = clt.read_ascii_column(cdir_in+'/cd_dtv_'+vtv_u[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])

xcd_s = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    vrt_s[jtv] = float(vtv_s[jtv])/100. # dTv as a float
    for ja in range(nb_algo):
        xcd_s[:,:,jtv,ja] = clt.read_ascii_column(cdir_in+'/cd_dtv_'+vtv_s[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])

xch_u = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    for ja in range(nb_algo):
        xch_u[:,:,jtv,ja] = clt.read_ascii_column(cdir_in+'/ch_dtv_'+vtv_u[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])

xch_s = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    for ja in range(nb_algo):
        xch_s[:,:,jtv,ja] = clt.read_ascii_column(cdir_in+'/ch_dtv_'+vtv_s[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])


xce_u = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    for ja in range(nb_algo):
        xce_u[:,:,jtv,ja] = clt.read_ascii_column(cdir_in+'/ce_dtv_'+vtv_u[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])

xce_s = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    for ja in range(nb_algo):
        xce_s[:,:,jtv,ja] = clt.read_ascii_column(cdir_in+'/ce_dtv_'+vtv_s[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])


##### L ####
xlo_u = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    #vrt_u[jtv] = float(vtv_u[jtv])/100. # dTv as a float
    for ja in range(nb_algo):
        xlo_u[:,:,jtv,ja] = clt.read_ascii_column(cdir_in+'/lo_dtv_'+vtv_u[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])

xlo_s = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    #vrt_s[jtv] = float(vtv_s[jtv])/100. # dTv as a float
    for ja in range(nb_algo):
        xlo_s[:,:,jtv,ja] = clt.read_ascii_column(cdir_in+'/lo_dtv_'+vtv_s[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])

##### UN10 ####
xun_u = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    #vrt_u[jtv] = float(vtv_u[jtv])/100. # dTv as a float
    for ja in range(nb_algo):
        xun_u[:,:,jtv,ja] = clt.read_ascii_column(cdir_in+'/un_dtv_'+vtv_u[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])
        xun_u[1,:,jtv,ja] = xun_u[1,:,jtv,ja] - xun_u[0,:,jtv,ja] ; # => UN10 - U

xun_s = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    #vrt_s[jtv] = float(vtv_s[jtv])/100. # dTv as a float
    for ja in range(nb_algo):
        xun_s[:,:,jtv,ja] = clt.read_ascii_column(cdir_in+'/un_dtv_'+vtv_s[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])
        xun_s[1,:,jtv,ja] = xun_s[1,:,jtv,ja] - xun_s[0,:,jtv,ja] ; # => UN10 - U

##### u* ####
xus_u = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    #vrt_u[jtv] = float(vtv_u[jtv])/100. # dTv as a float
    for ja in range(nb_algo):
        xus_u[:,:,jtv,ja] = clt.read_ascii_column(cdir_in+'/us_dtv_'+vtv_u[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])

xus_s = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    #vrt_s[jtv] = float(vtv_s[jtv])/100. # dTv as a float
    for ja in range(nb_algo):
        xus_s[:,:,jtv,ja] = clt.read_ascii_column(cdir_in+'/us_dtv_'+vtv_s[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])

##### Rib ####
xri_u = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    #vrt_u[jtv] = float(vtv_u[jtv])/100. # dTv as a float
    for ja in range(nb_algo):
        xri_u[:,:,jtv,ja] = clt.read_ascii_column(cdir_in+'/ri_dtv_'+vtv_u[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])

xri_s = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    #vrt_s[jtv] = float(vtv_s[jtv])/100. # dTv as a float
    for ja in range(nb_algo):
        xri_s[:,:,jtv,ja] = clt.read_ascii_column(cdir_in+'/ri_dtv_'+vtv_s[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])

##### z0 ####
xz0_u = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    #vrt_u[jtv] = float(vtv_u[jtv])/100. # dTv as a float
    for ja in range(nb_algo):
        xz0_u[:,:,jtv,ja] = clt.read_ascii_column(cdir_in+'/z0_dtv_'+vtv_u[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])

xz0_s = nmp.zeros((2,nU,ntv,nb_algo))
for jtv in range(ntv):
    #vrt_s[jtv] = float(vtv_s[jtv])/100. # dTv as a float
    for ja in range(nb_algo):
        xz0_s[:,:,jtv,ja] = clt.read_ascii_column(cdir_in+'/z0_dtv_'+vtv_s[jtv]+'_sst_'+csst+'_'+valgo_nm[ja]+'.dat', [0,1])




        

xy_labs = (-1.1, 1.)

dCd_u = 0.2 ; dCd_s = 0.5
dCe = 0.2

F_cd_u_min = 0.8
F_cd_u_max = 3.
F_cd_s_min = 0.
F_cd_s_max = 2.5

F_ch_u_min = 0.8
F_ch_u_max = 1.8
F_ch_s_min = 0.
F_ch_s_max = 1.4

F_ce_u_min = 0.8
F_ce_u_max = 1.8
F_ce_s_min = 0.
F_ce_s_max = 1.6



dF = 0.25

U_min = 0.
U_max = 34.
dU = 1.


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
def plot_Cx( nvt, nba, cvar_nm, vu10, xCx, vstab, valg_dn, w_min, w_max, dw, cx_min, cx_max, dcx, istab=0, cunit='' ):
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
    for jtv in range(nvt):
        for ja in range(nba):
            ax1.plot( vu10, xCx[1,:,jtv,ja], label=valg_dn[ja]+r' $\Delta\Theta=$'+csgn+str(vstab[jtv])+'K', \
                      linewidth=3.5-jtv*0.7, color=vcolor[ja], linestyle=vlines[ja], zorder=0.75 )
    vx_ticks = nmp.arange(w_min, w_max+dw, dw)
    plt.xticks( vx_ticks )
    locs, labels = plt.xticks() ; jl=0; newlabels = []
    for ll in locs:
        if not jl % 2 == 0: clab = ''
        else: clab = str(int(locs[jl]))
        newlabels.append(clab); jl=jl+1
    plt.xticks(locs,newlabels)
    ax1.set_xlim(U_min, U_max)
    vy_ticks = nmp.arange(cx_min, cx_max+dcx, dcx) ; plt.yticks( vy_ticks )
    ax1.set_ylim(cx_min, cx_max + dcx*0.8)
    #
    ctit = r'$'+cvar_nm+'$ --'+cstblt+'--'
    plt.text(0.5, y_tit, ctit, horizontalalignment='center', transform = ax1.transAxes, bbox=dict(boxstyle="square", fc='w'), **font_ttl)
    plt.ylabel(r''+cunit, **font_ylb)
    plt.xlabel(r'$U_{10m}$ [m/s]',     **font_xlb)
    plt.legend(bbox_to_anchor=(0.34, 0.16), ncol=3, shadow=False, fancybox=True)
    #
    plt.savefig(cvar_nm+'_'+cstblt+'.'+fig_frmt, dpi=DPI_FIG, facecolor='w', edgecolor='w', orientation='portrait')
    plt.close(1)
    return 0
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if not l_multi_fig:

    iplt = plot_Cx( ntv, nb_algo, 'C_D', vU, xcd_u, vrt_s, valgo_DN, U_min, U_max, dU, F_cd_u_min, F_cd_u_max, dCd_u, istab=0, cunit='$10^{-3}$' )
    iplt = plot_Cx( ntv, nb_algo, 'C_D', vU, xcd_s, vrt_s, valgo_DN, U_min, U_max, dU, F_cd_s_min, F_cd_s_max, dCd_u, istab=1, cunit='$10^{-3}$' )

    iplt = plot_Cx( ntv, nb_algo, 'C_H', vU, xch_u, vrt_s, valgo_DN, U_min, U_max, dU, F_ch_u_min, F_ch_u_max, dCd_u, istab=0, cunit='$10^{-3}$' )
    iplt = plot_Cx( ntv, nb_algo, 'C_H', vU, xch_s, vrt_s, valgo_DN, U_min, U_max, dU, F_ch_s_min, F_ch_s_max, dCd_u, istab=1, cunit='$10^{-3}$' )

    iplt = plot_Cx( ntv, nb_algo, 'C_E', vU, xce_u, vrt_s, valgo_DN, U_min, U_max, dU, F_ce_u_min, F_ce_u_max, dCd_u, istab=0, cunit='$10^{-3}$' )
    iplt = plot_Cx( ntv, nb_algo, 'C_E', vU, xce_s, vrt_s, valgo_DN, U_min, U_max, dU, F_ce_s_min, F_ce_s_max, dCd_u, istab=1, cunit='$10^{-3}$' )


    iplt = plot_Cx( ntv, nb_algo, 'L',   vU, xlo_u, vrt_s, valgo_DN, U_min, U_max, dU, -500.,   0., 50., istab=0, cunit='[m]' )
    iplt = plot_Cx( ntv, nb_algo, 'L',   vU, xlo_s, vrt_s, valgo_DN, U_min, U_max, dU,    0., 500., 50., istab=1, cunit='[m]' )

    iplt = plot_Cx( ntv, nb_algo, 'UN-U',   vU, xun_u, vrt_s, valgo_DN, U_min, U_max, dU,    0.,  2., 0.1, istab=0, cunit='[m/s]' )
    iplt = plot_Cx( ntv, nb_algo, 'UN-U',   vU, xun_s, vrt_s, valgo_DN, U_min, U_max, dU,   -2.,  0., 0.1, istab=1, cunit='[m/s]' )

    iplt = plot_Cx( ntv, nb_algo, 'u*',   vU, xus_u, vrt_s, valgo_DN, U_min, U_max, dU,    0.,  2., 0.1, istab=0, cunit='[m/s]' )
    iplt = plot_Cx( ntv, nb_algo, 'u*',   vU, xus_s, vrt_s, valgo_DN, U_min, U_max, dU,    0.,  2., 0.1, istab=1, cunit='[m/s]' )

    iplt = plot_Cx( ntv, nb_algo, 'Ri_B',   vU, xri_u, vrt_s, valgo_DN, U_min, U_max, dU,   -10.,   0., 1, istab=0, cunit='[]' )
    iplt = plot_Cx( ntv, nb_algo, 'Ri_B',   vU, xri_s, vrt_s, valgo_DN, U_min, U_max, dU,     0.,  10., 1, istab=1, cunit='[]' )

    iplt = plot_Cx( ntv, nb_algo, 'z_0',   vU, xz0_u, vrt_s, valgo_DN, U_min, U_max, dU,     0.,  0.1, 0.01, istab=0, cunit='[]' )
    iplt = plot_Cx( ntv, nb_algo, 'z_0',   vU, xz0_s, vrt_s, valgo_DN, U_min, U_max, dU,     0.,  0.1, 0.01, istab=1, cunit='[]' )

    
else:
    
    fig = plt.figure(num = 1, figsize=(12,14), facecolor='w', edgecolor='k')
    ax1  = plt.axes([rleft, 2*rheig+3*rbott+rbext, rwidt, rheig])
    ax1.grid(color='k', linestyle='-', linewidth=0.2, zorder=0.1)
    ax1.annotate('Unstable', xy=(0.4, 1.06), xycoords='axes fraction', **font_ttl)
    for jtv in range(ntv):
        for ja in range(nb_algo):
            ax1.plot(vU, xcd_u[1,:,jtv,ja], linewidth=3.5-jtv*0.7, color=vcolor[ja], zorder=0.75 )
    vx_ticks = nmp.arange(U_min, U_max+dU, dU)
    plt.xticks( vx_ticks )
    locs, labels = plt.xticks() ; jl=0; newlabels = []
    for ll in locs:
        if not jl % 2 == 0: clab = ''
        else: clab = str(int(locs[jl]))
        newlabels.append(clab); jl=jl+1
    plt.xticks(locs,newlabels)
    ax1.set_xlim(U_min, U_max)
    vy_ticks = nmp.arange(F_cd_u_min, F_cd_u_max+dCd_u, dCd_u) ; plt.yticks( vy_ticks )
    ax1.set_ylim(F_cd_u_min, F_cd_u_max + dCd_u*0.8)
    ax1.annotate('a)', xy=xy_labs, xycoords='axes fraction', **font_lab)
    ctit = r'$C_D$'
    plt.text(0.5, y_tit, ctit, horizontalalignment='center', transform = ax1.transAxes, bbox=dict(boxstyle="square", fc='w'), **font_ttl)
    if lshow_xylabs: plt.ylabel(r'$(10^{-3})$', **font_ylb)


    ########################################################################
    ax2  = plt.axes([0.5+rleft, 2*rheig+3*rbott+rbext, rwidt, rheig])
    ax2.grid(color='k', linestyle='-', linewidth=0.2, zorder=0.1)
    ax2.annotate('Stable', xy=(0.42, 1.06), xycoords='axes fraction', **font_ttl)
    for jtv in range(ntv):
        for ja in range(nb_algo):
            ax2.plot(vU, xcd_s[1,:,jtv,ja], label=valgo_DN[ja]+r'$\Delta\Theta=\pm$'+str(vrt_s[jtv])+'K', linewidth=3.5-jtv*0.7, color=vcolor[ja], zorder=0.75)
    vx_ticks = nmp.arange(U_min, U_max+dU, dU)
    plt.xticks( vx_ticks )
    plt.xticks(locs,newlabels)
    ax2.set_xlim(U_min, U_max)
    vy_ticks = nmp.arange(F_cd_s_min, F_cd_s_max+dCd_s, dCd_s) ; plt.yticks( vy_ticks )
    ax2.set_ylim(F_cd_s_min, F_cd_s_max + dCd_s/2)
    ax2.annotate('b)', xy=xy_labs, xycoords='axes fraction', **font_lab)
    ctit = r'$C_D$'
    plt.text(0.5, y_tit, ctit, horizontalalignment='center', transform = ax2.transAxes, bbox=dict(boxstyle="square", fc='w'), **font_ttl)



    ############## C_H ###################
    ax3  = plt.axes([rleft, rheig+2*rbott+rbext, rwidt, rheig])
    ax3.grid(color='k', linestyle='-', linewidth=0.2, zorder=0.1)
    for jtv in range(ntv):
        for ja in range(nb_algo):
            ax3.plot(vU, xch_u[1,:,jtv,ja], label=valgo_DN[ja]+r'$\Delta\Theta=$'+str(vrt_u[jtv]),  linewidth=3.5-jtv*0.7, color=vcolor[ja], zorder=0.75 )
    vx_ticks = nmp.arange(U_min, U_max+dU, dU)
    plt.xticks( vx_ticks )
    plt.xticks(locs,newlabels)
    ax3.set_xlim(U_min, U_max)
    vy_ticks = nmp.arange(F_ch_u_min, F_ce_u_max+dCe, dCe) ; plt.yticks( vy_ticks )
    ax3.set_ylim(F_ch_u_min, F_ch_u_max + dCe/2.)
    ax3.annotate('c)', xy=xy_labs, xycoords='axes fraction', **font_lab)

    ctit = r'$C_H$'
    plt.text(0.5, y_tit, ctit, horizontalalignment='center', transform = ax3.transAxes, bbox=dict(boxstyle="square", fc='w'), **font_ttl)
    if lshow_xylabs: plt.ylabel(r'$(10^{-3})$', **font_ylb)
    ax4  = plt.axes([0.5+rleft, rheig+2*rbott+rbext, rwidt, rheig])
    ax4.grid(color='k', linestyle='-', linewidth=0.2, zorder=0.1)
    for jtv in range(ntv):
        for ja in range(nb_algo):
            ax4.plot(vU, xch_s[1,:,jtv,ja], label=valgo_DN[ja]+r'$\Delta\Theta=$'+str(vrt_s[jtv]),  linewidth=3.5-jtv*0.7, color=vcolor[ja], zorder=0.75 )
    vx_ticks = nmp.arange(U_min, U_max+dU, dU)
    plt.xticks( vx_ticks )
    plt.xticks(locs,newlabels)
    ax4.set_xlim(U_min, U_max)
    vy_ticks = nmp.arange(F_ch_s_min, F_ce_s_max+dCe, dCe) ; plt.yticks( vy_ticks )
    ax4.set_ylim(F_ch_s_min, F_ch_s_max + dCe*0.8)
    ax4.annotate('d)', xy=xy_labs, xycoords='axes fraction', **font_lab)
    ctit = r'$C_H$'
    plt.text(0.5, y_tit, ctit, horizontalalignment='center', transform = ax4.transAxes, bbox=dict(boxstyle="square", fc='w'), **font_ttl)


    ############## C_E ###################
    ax5  = plt.axes([rleft, rbott+rbext, rwidt, rheig])
    ax5.grid(color='k', linestyle='-', linewidth=0.2, zorder=0.1)

    for jtv in range(ntv):
        for ja in range(nb_algo):
            ax5.plot(vU, xce_u[1,:,jtv,ja], label=valgo_DN[ja]+r'$\Delta\Theta=$'+str(vrt_u[jtv]),  linewidth=3.5-jtv*0.7, color=vcolor[ja], zorder=0.75 )

    vx_ticks = nmp.arange(U_min, U_max+dU, dU)
    plt.xticks( vx_ticks )
    plt.xticks(locs,newlabels)

    ax5.set_xlim(U_min, U_max)
    vy_ticks = nmp.arange(F_ce_u_min, F_ce_u_max+dCe, dCe) ; plt.yticks( vy_ticks )
    ax5.set_ylim(F_ce_u_min, F_ce_u_max + dCe/2.)

    ax5.annotate('e)', xy=xy_labs, xycoords='axes fraction', **font_lab)

    if lshow_xylabs: plt.xlabel(r'$U_{10m}$', **font_xlb)

    ctit = r'$C_E$'
    plt.text(0.5, y_tit, ctit, horizontalalignment='center', transform = ax5.transAxes, bbox=dict(boxstyle="square", fc='w'), **font_ttl)
    if lshow_xylabs: plt.ylabel(r'$(10^{-3})$', **font_ylb)
    ax6  = plt.axes([0.5+rleft, rbott+rbext, rwidt, rheig])

    ax6.grid(color='k', linestyle='-', linewidth=0.2, zorder=0.1)

    for jtv in range(ntv):
        for ja in range(nb_algo):
            ax6.plot(vU, xce_s[1,:,jtv,ja], label=valgo_DN[ja]+r'$\Delta\Theta=\pm$'+str(vrt_s[jtv])+'K', linewidth=3.5-jtv*0.7, color=vcolor[ja], zorder=0.75 )

    vx_ticks = nmp.arange(U_min, U_max+dU, dU)
    plt.xticks( vx_ticks )
    plt.xticks(locs,newlabels)

    ax6.set_xlim(U_min, U_max)
    vy_ticks = nmp.arange(F_ce_s_min, F_ce_s_max+dCe, dCe) ; plt.yticks( vy_ticks )
    ax6.set_ylim(F_ce_s_min, F_ce_s_max + dCe*0.8)
    ax6.annotate('f)', xy=xy_labs, xycoords='axes fraction', **font_lab)
    if lshow_xylabs: plt.xlabel(r'$U_{10m}$', **font_xlb)
    ctit = r'$C_E$'
    plt.text(0.5, y_tit, ctit, horizontalalignment='center', transform = ax6.transAxes, bbox=dict(boxstyle="square", fc='w'), **font_ttl)
    #plt.legend(loc='lower left', ncol=3, shadow=False, fancybox=True)
    plt.legend(bbox_to_anchor=(0.6, -0.22), ncol=3, shadow=False, fancybox=True)


    plt.savefig(cf_fig+'.'+fig_frmt, dpi=DPI_FIG, facecolor='w', edgecolor='w', orientation='portrait')
    plt.close(1)
