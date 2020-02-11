#!/bin/bash

YEAR=2018
DIR_IN="/home/datawork-lops-drakkarcom/DATA-REFERENCE/ERA5-FORCING/ROOT-FILES"

fo="ERA5_station_arctic_1h_${YEAR}.nc"

ip=1224
jp=18

# istl1
DIR_IN2="/mnt/meom/workdir/brodeau/ECMWF/ERA5-T720/ARCTIC/build"


if [ `hostname` = "ige-meom-cal1" ]; then

    for cv in "istl1" "skt"; do
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

for cv in "msl" "d2m" "t2m" "u10" "v10" "ssrd" "strd"; do
    echo
    if [ ! -f ./${cv}_era5.tmp ]; then
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
        CMD="ncks -4 -L 5 --cnk_dmn longitude,3 --cnk_dmn latitude,1 --cnk_dmn time,1 ${cv}_era5.tmp -o ${cv}_era5.nc4"
        echo; echo "${CMD}"; ${CMD}; echo
        rm -f ${cv}.tmp ${cv}0.tmp ${cv}_era5.tmp
        #
    else
        echo " ${cv}_era5.tmp is already here!"
    fi
    echo
done



if [ `hostname` = "ige-meom-cal1" ]; then
    rm -f ${fo}

    rsync -avP msl_era5.tmp ${fo}

    for cv in "d2m" "t2m" "u10" "v10" "ssrd" "strd"; do
        ncks -A --no-abc -h -v ${cv} ${cv}_era5.tmp -o ${fo}
    done
fi
