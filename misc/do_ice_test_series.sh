#!/bin/bash

if   [ "`hostname`" = "luitel" ]; then
    F_IN="/data/gcm_setup/STATION_ASF/sea-ice/ERA5_station_arctic_1h_2018_1x1.nc"
elif [ "`hostname`" = "merlat" ]; then
    F_IN="/MEDIA/data/STATION_ASF/ICE/ERA5_station_arctic_1h_2018_1x1.nc"
elif [ "`hostname`" = "salvelinus" ]; then
    F_IN="/home/laurent/data/ERA5_station_arctic_1h_2018_1x1.nc"
else
    echo "UNKNOW host!"; exit
fi

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

# Lupkes et al (2012)
./bin/test_aerobulk_ice_series.x -f ${F_IN} <<EOF
3
10
2
EOF

# Lupkes & Gryanik (2015)
./bin/test_aerobulk_ice_series.x -f ${F_IN} <<EOF
4
10
2
EOF


python3 ./python/plot_tests/plot_ice_bulk_comp.py .


