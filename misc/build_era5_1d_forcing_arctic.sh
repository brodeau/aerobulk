#!/bin/bash

YEAR=2018
DIR_IN="/home/datawork-lops-drakkarcom/DATA-REFERENCE/ERA5-FORCING/ROOT-FILES"

ip=1224
jp=18


file_tsk_ist=""

# ERA5_d2m_y2018.nc
# ERA5_msdwlwrf_y2018.nc
# ERA5_msdwswrf_y2018.nc
# ERA5_msl_y2018.nc
# ERA5_msr_y2018.nc
# ERA5_mtpr_y2018.nc
# ERA5_sf_y2018.nc
# ERA5_ssrd_y2018.nc
# ERA5_strd_y2018.nc
# ERA5_t2m_y2018.nc
# ERA5_tp_y2018.nc
# ERA5_u10_y2018.nc
# ERA5_v10_y2018.nc



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
