#!/bin/bash

#if   [ "`hostname`" = "luitel" ]; then
#    F_IN="/data/gcm_setup/STATION_ASF/sea-ice/ERA5_station_arctic_1h_2018_1x1.nc"
#elif [ "`hostname`" = "merlat" ]; then
#    F_IN="/MEDIA/data/AEROBULK/ERA5/Station_ERA5_81N_36p75E_ice_hourly_y2018.nc"; #REF!
#    #F_IN="/MEDIA/data/AEROBULK/ERA5/Station_ERA5_85p5N_54W_ice_hourly_y2018.nc"
#elif [ "`hostname`" = "ige-meom-cal1" ]; then
#    #F_IN="/mnt/meom/workdir/brodeau/AEROBULK/ERA5/Station_ERA5_81N_36p75E_ice_hourly_y2018.nc"; #REF!
#    F_IN="/mnt/meom/workdir/brodeau/AEROBULK/ERA5/Station_ERA5_85p5N_54W_ice_hourly_y2018.nc"
#elif [ "`hostname`" = "salvelinus" ]; then
#    F_IN="/home/laurent/data/ERA5_station_arctic_1h_2018_1x1.nc"
#else
#    echo "UNKNOW host!"; exit
#fi


F_IN="${HOME}/DEV/NEMO/NEMOGCM_trunk/tests/STATION_ASF/input_data/ERA5_NorthGreenland_surface_84N_-36E_1h_y2018.nc"


# NEMO default...
./bin/test_aerobulk_buoy_series_ice.x -f ${F_IN} -3 <<EOF
1
10
2
EOF

# Andreas et al (2005)
./bin/test_aerobulk_buoy_series_ice.x -f ${F_IN} -3 <<EOF
2
10
2
EOF

# Lupkes et al (2012)
./bin/test_aerobulk_buoy_series_ice.x -f ${F_IN} -3 <<EOF
3
10
2
EOF

# Lupkes & Gryanik (2015)
./bin/test_aerobulk_buoy_series_ice.x -f ${F_IN} -3 <<EOF
4
10
2
EOF


python3 ./python/plot_tests/plot_ice_bulk_comp.py .
