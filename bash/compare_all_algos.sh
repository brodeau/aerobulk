#!/bin/bash

export AEROBULK_HOME="${HOME}/DEV/aerobulk"

SST=15

cd ${AEROBULK_HOME}/

# First generating all the data with a SST of ${SST} C:
echo
echo " Will compute tests with 'test_cx_vs_wind.x' based on a SST of ${SST} deg.C!"

mkdir -p ${AEROBULK_HOME}/dat

for calgo in "coare3p6" "ncar" "ecmwf" "andreas"; do
    fdone="test_${calgo}_${SST}.done"
    if [ ! -f ${fdone} ]; then
        ./bin/test_cx_vs_wind.x ${calgo} ${SST}   1> out_test_${calgo}_${SST}.out 2> err_test_${calgo}_${SST}.err &
        touch ${fdone}
    else
        echo " *** ${calgo} has already been done..."
    fi
done
wait
echo

# Generating the plots:

./python/plot_tests/plot_Cx_wind.py

