CC      = gcc
CFLAGS  = -Wall -std=c99 -I../include -fopenmp -g

LDLIBS  = -lfftw3 -lgsl -lgslcblas -lm -fopenmp -g
OBJECTS = main.o wavelet.o CleanData.o ERSP.o convolution.o

morlet: $(OBJECTS)
	$(CC) -o $@  $(OBJECTS) $(LDLIBS)

clean:
	rm *.o morlet  *.log
