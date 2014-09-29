OBJS = hull.o
CFLAGS=`pkg-config --cflags --libs libgeodecomp`
CC=g++
DEBUG = -g

GOL_adcirc: $(OBJS) main.cpp 
	$(CC) $(CFLAGS) $(OBJS) main.cpp -o GOL_adc

clean: 
	rm *.o GOL_adc data/* plotting/plots/*