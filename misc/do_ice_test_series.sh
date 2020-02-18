#!/bin/bash

F_IN="/data/gcm_setup/STATION_ASF/sea-ice/ERA5_station_arctic_1h_2018_1x1.nc"

# NEMO default...
./bin/test_aerobulk_ice_series.x -f ${F_IN} <<EOF
1
10
2
EOF

# Andreas et al (2005)
./bin/test_aerobulk_ice_series.x -f ${F_IN} <<EOF
2
10
2
EOF

# Lupkes & Gryanik (2015)
./bin/test_aerobulk_ice_series.x -f ${F_IN} <<EOF
3
10
2
EOF
