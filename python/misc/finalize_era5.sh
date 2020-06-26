#!/bin/bash

fr=$1

ncrcat -O ${fr}m* -o ${fr}.nc

ncks -O --no-abc -4 -L 9 --cnk_dmn x,3 --cnk_dmn y,3 --cnk_dmn time_counter,365 ${fr}.nc -o ${fr}.nc4

mv -f ${fr}.nc ${fr}.nc.old

mv -f ${fr}.nc4 ${fr}.nc



