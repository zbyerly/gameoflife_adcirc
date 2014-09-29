#!/usr/bin/gnuplot
set terminal png size 1800,1800

unset key

#set xrange [-72.6 to -72.3]
#set yrange [40.8 to 40.9]

do for [t=0:100] {
set output sprintf('plots/output%i.png',t)
plot \
sprintf("../data/output0.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output1.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output2.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output3.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output4.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output5.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output6.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output7.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output8.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output9.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output10.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output11.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output12.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output13.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output14.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output15.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output16.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output17.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output18.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output19.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output20.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output21.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output22.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output23.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output24.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output25.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output26.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output27.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output28.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output29.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output30.%i.dat",t) u 2:3:4 w labels,\
sprintf("../data/output31.%i.dat",t) u 2:3:4 w labels
}
