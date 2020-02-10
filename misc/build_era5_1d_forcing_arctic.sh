#!/bin/bash

YEAR=2018
DIR_IN="/home/datawork-lops-drakkarcom/DATA-REFERENCE/ERA5-FORCING/ROOT-FILES"

fo="ERA5_station_arctic_1h_${YEAR}.nc"

ip=1224
jp=18

file_ist=""
file_tsk=""


for cv in "msl" "d2m" "t2m" "u10" "v10" "ssrd" "strd"; do
    echo
    if [ ! -f ./${cv}_era5.tmp ]; then
        #
        f_in=${DIR_IN}/${YEAR}/ERA5_${cv}_y${YEAR}.nc
        CMD="ncks -F -O --no-abc -h -d longitude,${ip},`expr ${ip} + 2` -d latitude,${jp},`expr ${jp} + 2` ${f_in} -o ${cv}.tmp"
        echo; echo "${CMD}"; ${CMD}; echo
        #
        ncecat -O ${cv}.tmp ${cv}_era5.tmp     # Add degenerate record dimension named "record"
        ncpdq -O -a time,record ${cv}_era5.tmp ${cv}_era5.tmp # Switch "record" and "time"
        ncwa -O -a record ${cv}_era5.tmp ${cv}_era5.tmp       # Remove (degenerate) "record"
        rm -f ${cv}.tmp
        #
    else
        echo " ${cv}_era5.tmp is already here!"
    fi
    echo
done




rm -f ${fo}

rsync -avP msl_era5.tmp ${fo}

for cv in "d2m" "t2m" "u10" "v10" "ssrd" "strd"; do
    ncks -A --no-abc -h -v ${cv} ${cv}_era5.tmp -o ${fo}
done

