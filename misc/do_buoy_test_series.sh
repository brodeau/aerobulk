#!/bin/bash

if [ `hostname` = "luitel"  ]; then DSTOR="/data/gcm_setup/STATION_ASF/STATION_ASF-I"; fi
if [ `hostname` = "lacroix" ]; then DSTOR="/data1/laurent/STATION_ASF/STATION_ASF-I" ; fi
if [ `hostname` = "merlat"  ]; then DSTOR="/MEDIA/data/STATION_ASF/STATION_ASF-I"    ; fi

for calgo in "ecmwf" "coare3p6"; do

    if   [ "${calgo}" = "ecmwf" ]; then
        id=2
    elif [ "${calgo}" = "coare3p6" ]; then
        id=3
    fi

    fout="out_test_aerobulk_buoy_${calgo}.out"
    
    rm -f lolo_${calgo}.nc ${fout}

    echo
    echo " *** Calling:"
    echo "  test_aerobulk_buoy_series_skin.x -f ${DSTOR}/Station_PAPA_50N-145W_2012-2012.nc4 1>${fout}"    
    ./bin/test_aerobulk_buoy_series_skin.x -f ${DSTOR}/Station_PAPA_50N-145W_2012-2012.nc4 1>${fout} <<EOF
${id}
14
14
EOF

    ncks -O -d time,3500,4500 lolo_${calgo}.nc -o fluxes_short_${calgo}.nc

    echo "Done!"
    echo
    
done


