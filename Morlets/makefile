CC        = gcc
CXX       = g++
CFLAGS    = -fopenmp
CPPFLAGS  = -Wall -I../include -g -DNO_FREETYPE

LDLIBS  = -lfftw3 -lgsl -lgslcblas -lm -fopenmp -lpngwriter -lz
OBJECTS = main.o wavelet.o CleanData.o ERSP.o convolution.o Plot.o

morlet: $(OBJECTS)
	$(CXX) -o $@  $(OBJECTS) $(LDLIBS)

clean:
	rm *.o morlet  *.log
