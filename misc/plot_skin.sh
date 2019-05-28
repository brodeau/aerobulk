#!/bin/bash

RAD_LW=300

# Night unstable:
###################

RAD_SW=0
SST=22
Tair=20

# COARE:
./bin/test_skin_corr.x <<EOF
1
${RAD_SW}
${RAD_LW}
22
${Tair}
EOF

F1="dat/dT_skin_vs_wind_coare_SST${SST}_SW`printf "%04d" ${RAD_SW}`_LW`printf "%04d" ${RAD_LW}`_RH80_unstable.dat"

# ECMWF:
./bin/test_skin_corr.x <<EOF
2
${RAD_SW}
${RAD_LW}
${SST}
${Tair}
EOF

F2="dat/dT_skin_vs_wind_ecmwf_SST${SST}_SW`printf "%04d" ${RAD_SW}`_LW`printf "%04d" ${RAD_LW}`_RH80_unstable.dat"


FEPS="comp_SST${SST}_SW`printf "%04d" ${RAD_SW}`_LW`printf "%04d" ${RAD_LW}`_RH80_unstable.eps"

cat > plot.gp <<EOF
set size 1. , 1.6
set output '${FEPS}'
#set terminal postscript eps color portrait enhanced 'OpenSans' 18
set terminal postscript eps color enhanced 'OpenSans' 18


set ylabel 'T_{skin} - SST [K]'
set xlabel 'Wind speed [m/s]'
#set title  'depth_t(jpk)' font 'arial bold,20'

set ytic 0.1
set grid
#set key bottom left Left reverse
set key bottom left Left reverse box


set xrange [0:30]
set yrange [-1.:1.]

plot 0 notitle w l lt -1 lw 2, \\
     '${F1}' t 'COARE 3.0' w l lt 1 lw 3, \\
     '${F2}' t 'ECMWF'     w l lt 3 lw 3
EOF
gnuplot plot.gp
exit


#     '${F1}' t 'COARE 3.0' w lp lt 1 lw 2 ps 1.5 pt 6, \\
#     '${F2}' t 'ECMWF'     w lp lt 3 lw 2 ps 1.5 pt 8










# Noon unstable:
################


RAD_SW=900

# COARE:
./bin/test_skin_corr.x <<EOF
1
${RAD_SW}
${RAD_LW}
${SST}
${Tair}
EOF



# ECMWF:
./bin/test_skin_corr.x <<EOF
2
${RAD_SW}
${RAD_LW}
${SST}
${Tair}
EOF








