#!/bin/bash

# coare3p6

fout="out_test_aerobulk_buoy_coare3p6.out"


DSTOR="/data/gcm_setup/STATION_ASF/STATION_ASF-I"

if [ `hostname` = "lacroix" ]; then DSTOR="/data1/laurent/STATION_ASF/STATION_ASF-I" ; fi

rm -f lolo_coare3p6.nc ${fout}

./bin/test_aerobulk_buoy_series_skin.x -f ${DSTOR}/Station_PAPA_50N-145W_2012-2012.nc4 1>${fout} <<EOF
3
14
14
EOF


