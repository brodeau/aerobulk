#!/bin/bash

# ecmwfn

fout="out_test_aerobulk_buoy_ecmwfn.out"

if [ `hostname` = "luitel"  ]; then DSTOR="/data/gcm_setup/STATION_ASF/STATION_ASF-I"; fi
if [ `hostname` = "lacroix" ]; then DSTOR="/data1/laurent/STATION_ASF/STATION_ASF-I" ; fi
if [ `hostname` = "merlat"  ]; then DSTOR="/MEDIA/data/STATION_ASF/STATION_ASF-I"    ; fi

rm -f lolo_ecmwfn.nc ${fout}

./bin/test_aerobulk_buoy_series_skin.x -f ${DSTOR}/Station_PAPA_50N-145W_2012-2012_short.nc4 1>${fout} <<EOF
5
14
14
EOF

ncks -O -d time,3500,4500 lolo_ecmwfn.nc -o short_e.nc
