#!/bin/bash


# Air temperature of -3 deg.C:
./bin/test_ice.x <<EOF
-3
EOF

FOUT="z0_z0t_z0q__ustar_test.dat"











FEPS=`echo ${FOUT} | sed -e "s|.dat|.eps|g"`

cat > plot.gp <<EOF
set size 2. , 1.6
set output '${FEPS}'
#set terminal postscript eps color portrait enhanced 'OpenSans' 18
set terminal postscript eps color enhanced 'OpenSans' 18


set logscale y 10

set ylabel 'z_0 [m]'
set xlabel 'u* [m/s]'
#set title  'depth_t(jpk)' font 'arial bold,20'

set grid
#set key bottom left Left reverse
set key top right Left reverse box


set xrange [0:0.7]

set yrange [1.E-6:1.E-1]
set format y "10^{%L}"
#set ytic 0.0001

plot '${FOUT}' u 1:2 t ' z_0  ' w l lt rgb "#000000" lw 6, \\
     '${FOUT}' u 1:3 t ' z_0t ' w l lt rgb "#008ab8" lw 6, \\
     '${FOUT}' u 1:4 t ' z_0q ' w l lt rgb "#ffed00" lw 6
EOF
gnuplot plot.gp

#     '${F2}' t 'ECMWF'     w l lt rgb "#008ab8" lw 6, \\
#     '${F3}' notitle       w l lt rgb "#ffed00" lw 5, \\
#     '${F4}' notitle       w l lt rgb "#008ab8" lw 5



















list=`\ls *.eps`
for ff in ${list}; do
    fn=`echo ${ff} | sed -e s/".eps"/".png"/g`
    CMD="convert -density 180 -flatten ${ff} ${fn}"
    echo ${CMD}
    ${CMD}
done







exit


#     '${F1}' t 'COARE 3.0' w lp lt 1 lw 2 ps 1.5 pt 6, \\
#     '${F2}' t 'ECMWF'     w lp lt 3 lw 2 ps 1.5 pt 8
