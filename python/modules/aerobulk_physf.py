# Misc :

import sys
import numpy as nmp
import math

rt0 = 273.16


grav  = 9.8          # gravity
Rgas  = 287.04     
Patm  = 101000.    
ctv   = 0.608        # for virtual temperature


R_dry = 287.05       # Specific gas constant for dry air              [J/K/kg]
R_vap = 461.495      # Specific gas constant for water vapor          [J/K/kg]
reps0 = R_dry/R_vap  # ratio of gas constant for dry air and water vapor => ~ 0.622

cte   = 0.622     
kappa = 0.4          # Von Karman's constant
Cp    = 1000.5    
Pi    = 3.141592654 
eps_w = 0.987        # emissivity of water
sigma = 5.67E-8      # Stefan Boltzman constamt
alfa  = 0.066        # Surface albedo over ocean

rtt0 = 273.16     # triple point of temperature    [K]

sensit = 0.1

def Lvap(zsst):
    #
    # INPUT  : zsst => water temperature in [K]
    # OUTPUT : Lvap => Latent Heat of Vaporization [J/Kg]
    return ( 2.501 - 0.00237*(zsst - rt0) )*1.E6

def e_air(q_air, zslp):
    #
    #--------------------------------------------------------------------
    #                  **** Function e_air ****
    #
    # Gives vapour pressure of air from pressure and specific humidity
    #
    #--------------------------------------------------------------------
    #

    diff  = 1.E8
    e_old = q_air*zslp/reps0

    while diff > 1.:
        #print "Again... diff = ", diff
        ee = q_air/reps0*(zslp - (1. - reps0)*e_old)
        diff  = nmp.sum(abs( ee - e_old ))
        e_old = ee

    return ee





### Update: June 2019, LB:

def e_sat(rT):
    #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # rT:     air temperature [K]
    # e_sat:  water vapor at saturation [Pa]
    #
    # Recommended by WMO
    #
    # Goff, J. A., 1957: Saturation pressure of water on the new kelvin
    # temperature scale. Transactions of the American society of heating
    # and ventilating engineers, 347--354.
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #
    ztmp = nmp.zeros(nmp.shape(rT))
    ztmp = rtt0/rT
    #
    #e_sat = 100.*( 10.^(10.79574*(1. - ztmp) - 5.028*LOG10(rT/rtt0) \
    #                    + 1.50475*10.^(-4)*(1. - 10.^(-8.2969*(rT/rtt0 - 1.)) ) \
    #                    + 0.42873*10.^(-3)*(10.^(4.76955*(1. - ztmp)) - 1.) + 0.78614) )
    #
    e_sat = 100.*( nmp.power(10.,(10.79574*(1. - ztmp) - 5.028*nmp.log10(1./ztmp) \
                   + 1.50475*0.0001*(1. - nmp.power(10.,(-8.2969*(1./ztmp - 1.))) ) \
                   + 0.42873*0.001 *(nmp.power(10.,(4.76955*(1. - ztmp))) - 1.) + 0.78614 ) ) )
    #
    del ztmp
    #
    return e_sat

def q_air_dp(da, slp):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Air specific humidity from dew point temperature
    #     da          !: dew-point temperature   [K]
    #     slp         !: atmospheric pressure    [Pa]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #
    es = e_sat(da)
    q_air_dp = es*reps0/(slp - (1. - reps0)*es)
    return q_air_dp
