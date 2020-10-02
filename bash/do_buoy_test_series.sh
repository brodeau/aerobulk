#!/bin/bash

if [ `hostname` = "luitel"  ]; then DSTOR="/data/gcm_setup/STATION_ASF/aerobulk"; fi
if [ `hostname` = "lacroix" ]; then DSTOR="/data1/laurent/STATION_ASF/aerobulk" ; fi
if [ `hostname` = "merlat"  ]; then DSTOR="/MEDIA/data/STATION_ASF/aerobulk"    ; fi


#fforcing="Station_PAPA_50N-145W_2012-2012.nc4" ; zu=14 ; zt=14      ; clabel="PAPA_50N-145W"
fforcing="idealized_forcing_test_STATION_ASF_1h.nc4" ; zu=10 ; zt=2  ; clabel="IDEALIZED"


for calgo in "ecmwf" "coare3p6"; do

    if   [ "${calgo}" = "ecmwf" ]; then
        id=2
    elif [ "${calgo}" = "coare3p6" ]; then
        id=3
    fi

    fout="out_test_aerobulk_buoy_${calgo}.out"
    
    rm -f ${clabel}_${calgo}.nc ${fout}

    echo
    echo " *** Calling:"
    echo "  test_aerobulk_buoy_series_skin.x -f ${DSTOR}/${fforcing} 1>${fout}"    
    ./bin/test_aerobulk_buoy_series_skin.x -f ${DSTOR}/${fforcing} -n ${clabel}  1>${fout} <<EOF
${id}
${zu}
${zt}
EOF

    #ncks -O -d time,3500,4500 ${clabel}_${calgo}.nc -o fluxes_short_${calgo}.nc

    echo "Done!"
    echo
    
done


