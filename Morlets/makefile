CC       = gcc
CXX      = g++
CFLAGS   = -I/opt/local/lib/gcc5/gcc/x86_64-apple-darwin15/5.4.0/include
CXXFLAGS = -Wall -Wextra -g -fopenmp -pedantic-errors `freetype-config --cflags` -L/opt/local/lib/

LDLIBS  = -lfftw3 -lgsl -lgslcblas -lm -lz `freetype-config --libs` -lz -lpngwriter -lpng
OBJECTS = main.o wavelet.o CleanData.o ERSP.o convolution.o Plot.o

morlet: $(OBJECTS)
	$(CXX) -o $@  $(OBJECTS) $(LDLIBS);
	install_name_tool -change libpngwriter.dylib /opt/local/lib/libpngwriter.dylib morlet

clean:
	rm *.o morlet  *.log
