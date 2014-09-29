#!/usr/bin/gnuplot
set terminal png

#set xrange [20000 to 120000]
#set yrange [-10000 to 120000]

do for [t=0:50] {
set output sprintf('output%i.png',t)
plot \
sprintf("output0.%i.dat",t) u 2:3:4 w labels,\
sprintf("output1.%i.dat",t) u 2:3:4 w labels,\
sprintf("output2.%i.dat",t) u 2:3:4 w labels,\
sprintf("output3.%i.dat",t) u 2:3:4 w labels
}