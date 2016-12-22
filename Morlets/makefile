CC      = g++
CFLAGS  = -Wall -I../include -fopenmp -g -lpng -lpngwriter -lz -DNO_FREETYPE

LDLIBS  = -lfftw3 -lgsl -lgslcblas -lm -fopenmp -L /Users/vinay/usr/local/lib
OBJECTS = main.o wavelet.o CleanData.o ERSP.o convolution.o Plot.o

morlet: $(OBJECTS)
	$(CC) -o $@  $(OBJECTS) $(LDLIBS)

clean:
	rm *.o morlet  *.log
