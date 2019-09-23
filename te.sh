#!/bin/bash

# ecmwf

fout="out_test_aerobulk_buoy_ecmwf.out"

if [ `hostname` = "luitel"  ]; then DSTOR="/data/gcm_setup/STATION_ASF/STATION_ASF-I"; fi
if [ `hostname` = "lacroix" ]; then DSTOR="/data1/laurent/STATION_ASF/STATION_ASF-I" ; fi
if [ `hostname` = "merlat"  ]; then DSTOR="/MEDIA/data/STATION_ASF/STATION_ASF-I"    ; fi

rm -f lolo_ecmwf.nc ${fout}

./bin/test_aerobulk_buoy_series_skin.x -f ${DSTOR}/Station_PAPA_50N-145W_2012-2012.nc4 1>${fout} <<EOF
2
14
14
EOF


