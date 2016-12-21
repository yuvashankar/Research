CC      = g++
CFLAGS  = -Wall -std=c99 -I../include -fopenmp -g -lpng -lpngwriter -lz -DNO_FREETYPE

LDLIBS  = -lfftw3 -lgsl -lgslcblas -lm -fopenmp -L../Users/vinay/build/lib/
OBJECTS = main.o wavelet.o CleanData.o ERSP.o convolution.o Plot.o

morlet: $(OBJECTS)
	$(CC) -o $@  $(OBJECTS) $(LDLIBS)

clean:
	rm *.o morlet  *.log
