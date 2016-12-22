#CC        = gcc
CXX       = gcc-5
#CFLAGS    = -fopenmp
CXXFLAGS  = -Wall -g -fopenmp `freetype-config --cflags` -std=c99

LDLIBS  = -lfftw3 -lgsl -lgslcblas -lm -lz `freetype-config --libs` -lz -lpngwriter -lpng -L../lib
OBJECTS = main.o wavelet.o CleanData.o ERSP.o convolution.o Plot.o

morlet: $(OBJECTS)
	$(CXX) -o $@  $(OBJECTS) $(LDLIBS)

clean:
	rm *.o morlet  *.log
