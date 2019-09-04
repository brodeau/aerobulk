#!/bin/bash

rm -f lolo.nc out_test_aerobulk_buoy.out

./bin/test_aerobulk_buoy_series_skin.x -f /data/gcm_setup/STATION_ASF/STATION_ASF-I/Station_PAPA_50N-145W_2012-2012.nc4 1>out_test_aerobulk_buoy.out <<EOF
14
14
EOF


