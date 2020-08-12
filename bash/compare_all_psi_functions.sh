#!/bin/bash

export AEROBULK_HOME="${HOME}/DEV/aerobulk"

SST=15

cd ${AEROBULK_HOME}/

# First generating the data:
echo
echo " Generating the data with 'test_psi_stab.x' !"
./bin/test_psi_stab.x
echo

# Generating the plots:
echo
echo " Generating the plots with 'plot_Psi_profiles.py' !"
./python/plot_tests/plot_Psi_profiles.py psi.nc




