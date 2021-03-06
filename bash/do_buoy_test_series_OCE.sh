#!/bin/bash

if [ `hostname` = "luitel"  ]; then DSTOR="/data/gcm_setup/STATION_ASF/aerobulk"; fi
if [ `hostname` = "lacroix" ]; then DSTOR="/data1/laurent/STATION_ASF/aerobulk" ; fi
if [ `hostname` = "merlat"  ]; then DSTOR="/MEDIA/data/STATION_ASF/aerobulk"    ; fi


fforcing="Station_PAPA_50N-145W_2012-2012.nc4" ; zu=14 ; zt=14      ; clabel="PAPA_50N-145W"
#fforcing="idealized_forcing_test_STATION_ASF_1h.nc4" ; zu=10 ; zt=2  ; clabel="IDEALIZED"


for calgo in "ecmwf" "coare3p6" "coare3p0" "ncar" "andreas" ; do

# "ncar" => 0 , "coare3p0" => 1 , "ecmwf" => 2 , "coare3p6" => 3, "andreas" => 4
    
    if   [ "${calgo}" = "ncar" ]; then
        id=0
    elif [ "${calgo}" = "coare3p0" ]; then
        id=1
    elif [ "${calgo}" = "ecmwf" ]; then
        id=2
    elif [ "${calgo}" = "coare3p6" ]; then
        id=3
    elif [ "${calgo}" = "andreas" ]; then
        id=4
    else
        echo "UNKNOWN ALGO: ${calgo} !"; exit
    fi

    fout="out_test_aerobulk_buoy_${calgo}.out"
    
    rm -f ${clabel}_${calgo}.nc ${fout}
    
    echo
    CMD="./bin/test_aerobulk_buoy_series_oce.x -f ${DSTOR}/${fforcing} -n ${clabel} -r -w"
    echo " *** Calling:"
    echo "  ${CMD}"
    ${CMD} 1>${fout} <<EOF
${id}
${zu}
${zt}
EOF

    #ncks -O -d time,3500,4500 ${clabel}_${calgo}.nc -o fluxes_short_${calgo}.nc

    echo "Done!"
    echo
    
done


