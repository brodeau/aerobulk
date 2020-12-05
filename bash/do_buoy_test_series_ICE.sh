#!/bin/bash

#if [ `hostname` = "luitel"  ]; then DSTOR="/data/gcm_setup/STATION_ASF/aerobulk"; fi
#if [ `hostname` = "lacroix" ]; then DSTOR="/data1/laurent/STATION_ASF/aerobulk" ; fi
#if [ `hostname` = "merlat"  ]; then DSTOR="/MEDIA/data/STATION_ASF/aerobulk"    ; fi

#


#fforcing="${HOME}/DEV/NEMO/NEMOGCM_trunk/tests/STATION_ASF/input_data/ERA5_arctic_surface_81N_36p75E_1h_y2018.nc" ; FORCING="ERA5_arctic"
fforcing="${HOME}/DEV/NEMO/NEMOGCM_trunk/tests/STATION_ASF/input_data/ERA5_NorthGreenland_surface_84N_-36E_1h_y2018.nc" ; FORCING="ERA5_NorthGreenland"
zu="10"
zt="2"

#  * "NEMO default (v4)"       => 1
#  * "Andreas (2005)"          => 2
#  * "Lupkes et al (2012)"     => 3
#  * "Lupkes & Gryanik (2015)" => 4


for calgo in "nemo" "an05" "lu12" "lg15"; do
    
    if   [ "${calgo}" = "nemo" ]; then
        id=1
    elif [ "${calgo}" = "an05" ]; then
        id=2
    elif [ "${calgo}" = "lu12" ]; then
        id=3
    elif [ "${calgo}" = "lg15" ]; then
        id=4
    else
        echo "UNKNOWN ALGO: ${calgo} !"; exit
    fi

    fout="out_test_aerobulk_buoy_${calgo}.out"
    
    rm -f ${clabel}_${calgo}.nc ${fout}
    
    echo
    CMD="./bin/test_aerobulk_buoy_series_ice.x -f ${fforcing} -n ${FORCING} -3"
    echo " *** Calling:"
    echo "  ${CMD}"
    ${CMD} 1>${fout} <<EOF
${id}
${zu}
${zt}
EOF

    echo "Done!"
    echo

done
