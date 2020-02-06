#!/bin/bash

./bin/test_aerobulk_buoy_series_ice.x -f /data/gcm_setup/STATION_ASF/sea-ice/NGreenLand_ERA5_Arctic_201901_1h_1D.nc4 <<EOF
0
10
2
EOF


./bin/test_aerobulk_buoy_series_ice.x -f /data/gcm_setup/STATION_ASF/sea-ice/NGreenLand_ERA5_Arctic_201901_1h_1D.nc4 <<EOF
1
10
2
EOF


