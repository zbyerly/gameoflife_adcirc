#!/usr/bin/gnuplot
set terminal png
set output "hulltest.png"
plot \
'points0.dat' w p, \
'points1.dat' w p, \
'points2.dat' w p, \
'points3.dat' w p, \
'hull0.dat' w l, \
'hull1.dat' w l, \
'hull2.dat' w l, \
'hull3.dat' w l
