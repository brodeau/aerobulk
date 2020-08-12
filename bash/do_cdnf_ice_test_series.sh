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

# LU12 general
./bin/test_aerobulk_cdnf_series.x -f ${F_IN} <<EOF
1
10
EOF

# LU12 light
./bin/test_aerobulk_cdnf_series.x -f ${F_IN} <<EOF
2
10
EOF

# LU13
./bin/test_aerobulk_cdnf_series.x -f ${F_IN} <<EOF
3
10
EOF

# LU15
./bin/test_aerobulk_cdnf_series.x -f ${F_IN} <<EOF
4
10
EOF

# LU15 light
./bin/test_aerobulk_cdnf_series.x -f ${F_IN} <<EOF
5
10
EOF
