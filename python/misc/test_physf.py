#!/usr/bin/env python
#
# L. Brodeau, Feb.2001

import sys
import numpy as nmp
#from netCDF4 import Dataset
#import string
#from os.path import exists,basename

import aerobulk_physf as abp


rd2m = nmp.zeros((2,1))
rslp = nmp.zeros((2,1))


rd2m[0,0] = 293.15 # 20 C
rslp[0,0] = 101000.0

rd2m[1,0] = 280.0
rslp[1,0] = 101000.0

print ' *** Dew-point temp. =', rd2m
print ' *** SLP =',             rslp


rq2m = abp.q_air_dp(rd2m, rslp)               


print '\n  => q2m =', 1000. * rq2m, 'g/kg\n'
