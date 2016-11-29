CC      = gcc
CFLAGS  = -Wall -O2 -std=c99 -I../include -fopenmp

LDLIBS  = -lfftw3 -lgsl -lgslcblas -lm -fopenmp
OBJECTS = main.o wavelet.o CleanData.o ERSP.o convolution.o

morlet: $(OBJECTS)
	$(CC) -o $@  $(OBJECTS) $(LDLIBS)

clean:
	rm *.o morlet  *.log
