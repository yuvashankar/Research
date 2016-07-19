CC      = gcc
CFLAGS  = -Wall -g -O2 -fopenmp -std=c99 -I../../bin/include -I/usr/include/openmpi-x86_64 -I/usr/local/cuda/include

LDFLAGS = -L../../bin/lib
LDLIBS  = -lfftw3 -lgsl -lgslcblas -lm -fopenmp
OBJECTS = main.o wavelet.o

morlet: $(OBJECTS)
	$(CC) -o $@  $(OBJECTS) $(LDFLAGS) $(LDLIBS)

clean:
	rm *.o morlet  *.log
