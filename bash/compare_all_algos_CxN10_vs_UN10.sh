#!/bin/bash

export AEROBULK_HOME="${HOME}/DEV/aerobulk"

cd ${AEROBULK_HOME}/

# First generating all the data with a SST of ${SST} C:
echo
echo " Will compute tests with 'bin/test_coef_n10.x'!"

mkdir -p ${AEROBULK_HOME}/dat

echo

./bin/test_coef_n10.x

echo

#exit

# Generating the plots:

python3 ./python/plot_tests/plot_CxN10_UN10.py

