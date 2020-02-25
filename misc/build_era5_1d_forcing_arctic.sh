#!/bin/bash

YEAR=2018
DIR_IN="/home/datawork-lops-drakkarcom/DATA-REFERENCE/ERA5-FORCING/ROOT-FILES"

fo="ERA5_station_arctic_1h_${YEAR}.nc"

coord="85.5 N, -54.0 W"
ip=1224 ; clon="306.0"
jp=18   ; clat="85.5"

# istl1
DIR_IN2="/mnt/meom/workdir/brodeau/ECMWF/ERA5-T720/ARCTIC/build"


if [ `hostname` = "ige-meom-cal1" ]; then

    #for cv in "istl1" "skt"; do
    for cv in "siconc"; do
        echo
        ftmp="./${cv}_era5.nc4"
        #
        if [ ! -f ${ftmp} ]; then
            #
            f_in=${DIR_IN2}/ERA5_${cv}_y${YEAR}.nc
            CMD="ncks -F -O --no-abc -h -d longitude,${ip},`expr ${ip} + 2` -d latitude,${jp},`expr ${jp} + 2` ${f_in} -o ${cv}0.tmp"
            echo; echo "${CMD}"; ${CMD}; echo
            #
            CMD="ncpdq -O -U ${cv}0.tmp -o ${cv}.tmp"
            echo; echo "${CMD}"; ${CMD}; echo
            #
            ncecat -O ${cv}.tmp ${cv}_era5.tmp     # Add degenerate record dimension named "record"
            ncpdq -O -a time,record ${cv}_era5.tmp ${cv}_era5.tmp # Switch "record" and "time"
            ncwa -O -a record ${cv}_era5.tmp ${cv}_era5.tmp       # Remove (degenerate) "record"
            #
            CMD="ncks --no-abc -4 -L 5 --cnk_dmn longitude,3 --cnk_dmn latitude,1 --cnk_dmn time,1 ${cv}_era5.tmp -o ${ftmp}"
            echo; echo "${CMD}"; ${CMD}; echo
            #
            #rm -f ${cv}.tmp ${cv}0.tmp ${cv}_era5.tmp
            #
        else
            echo " ${cv}_era5.tmp is already here!"
        fi
        echo
    done

fi


if [ `hostname` != "ige-meom-cal1" ]; then

    for cv in "msl" "d2m" "t2m" "u10" "v10" "ssrd" "strd"; do
        #
        ftmp="./${cv}_era5.nc4"
        #
        if [ ! -f ${ftmp} ]; then
            #
            f_in=${DIR_IN}/${YEAR}/ERA5_${cv}_y${YEAR}.nc
            CMD="ncks -F -O --no-abc -h -d longitude,${ip},`expr ${ip} + 2` -d latitude,${jp},`expr ${jp} + 2` ${f_in} -o ${cv}0.tmp"
            echo; echo "${CMD}"; ${CMD}; echo
            #
            CMD="ncpdq -O -U ${cv}0.tmp -o ${cv}.tmp"
            echo; echo "${CMD}"; ${CMD}; echo
            #
            ncecat -O ${cv}.tmp ${cv}_era5.tmp     # Add degenerate record dimension named "record"
            ncpdq -O -a time,record ${cv}_era5.tmp ${cv}_era5.tmp # Switch "record" and "time"
            ncwa -O -a record ${cv}_era5.tmp ${cv}_era5.tmp       # Remove (degenerate) "record"
            #
            #
            if [ "${cv}" = "ssrd" ] || [ "${cv}" = "strd" ]; then
                mv -f ${cv}_era5.tmp aaa.tmp
                ncap2 -h -O -s "${cv}0=${cv}/3600." aaa.tmp -o bbb.tmp
                ncks -h -v "${cv}0" bbb.tmp -o ${cv}_era5.tmp
                ncrename -v ${cv}0,${cv} ${cv}_era5.tmp
                ncatted -h -O -a units,${cv},o,c,"W m**-2" ${cv}_era5.tmp
                rm -f aaa.tmp bbb.tmp
            fi
            #
            #
            CMD="ncks -4 -L 5 --cnk_dmn longitude,3 --cnk_dmn latitude,1 --cnk_dmn time,1 ${cv}_era5.tmp -o ${ftmp}"
            echo; echo "${CMD}"; ${CMD}; echo
            rm -f ${cv}.tmp ${cv}0.tmp ${cv}_era5.tmp
            #
        else
            echo " ${cv}_era5.tmp is already here!"
        fi
        echo
    done

fi


exit

if [ `hostname` = "ige-meom-cal1" ]; then
    rm -f ${fo}

    rsync -avP msl_era5.nc4 ${fo}

    for cv in "istl1" "skt" "d2m" "t2m" "u10" "v10" "ssrd" "strd"; do
        ncks -A --no-abc -h -v ${cv} ${cv}_era5.nc4 -o ${fo}
    done

    ncrename -h -v ssrd,radsw ${fo}
    ncrename -h -v strd,radlw ${fo}
    
    ncatted -h -O -a About,global,o,c,"Extraction of ERA5 at ${coord} (i=${ip}, j=${jp})" ${fo}

fi
