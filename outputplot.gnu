#!/usr/bin/gnuplot
set terminal png
set output "output.png"

set xrange [20000 to 120000]
set yrange [-10000 to 120000]

plot \
"output0.0.dat" u 2:3:1 w labels