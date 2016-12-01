CC      = gcc
CFLAGS  = -Wall -std=c99 -I../include -fopenmp -g -DN_DEBUG

LDLIBS  = -lfftw3 -lgsl -lgslcblas -lm -fopenmp
OBJECTS = main.o wavelet.o CleanData.o ERSP.o convolution.o

morlet: $(OBJECTS)
	$(CC) -o $@  $(OBJECTS) $(LDLIBS)

clean:
	rm *.o morlet  *.log
