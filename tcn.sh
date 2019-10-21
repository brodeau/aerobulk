#!/bin/bash

# coare3p6n

fout="out_test_aerobulk_buoy_coare3p6n.out"

if [ `hostname` = "luitel"  ]; then DSTOR="/data/gcm_setup/STATION_ASF/STATION_ASF-I"; fi
if [ `hostname` = "lacroix" ]; then DSTOR="/data1/laurent/STATION_ASF/STATION_ASF-I" ; fi
if [ `hostname` = "merlat"  ]; then DSTOR="/MEDIA/data/STATION_ASF/STATION_ASF-I"    ; fi

rm -f lolo_coare3p6.nc ${fout} short.nc aa.out


###./bin/test_aerobulk_buoy_series_skin.x -f ${DSTOR}/Station_PAPA_50N-145W_2012-2012.nc4 1>${fout} <<EOF

./bin/test_aerobulk_buoy_series_skin.x -f ${DSTOR}/Station_PAPA_50N-145W_2012-2012_short.nc4 1>${fout} <<EOF
4
14
14
EOF


cat ${fout} | grep '#LBD' > aa.out
#ncks -O -d time,3500,4500 lolo_coare3p6.nc -o short.nc
