CC       = gcc
CXX      = g++
CFLAGS   = -I/opt/local/lib/gcc5/gcc/x86_64-apple-darwin15/5.4.0/include -I../include -Wall -Wextra -g -fopenmp -pedantic-errors -DNDEBUG
CXXFLAGS = -Wall -Wextra -g -fopenmp `freetype-config --cflags` -I/opt/local/lib/gcc5/gcc/x86_64-apple-darwin15/5.4.0/include -I../include -DNDEBUG

LDLIBS  = -lfftw3 -lfftw3_threads -lgsl -lgslcblas -lm -lz -lpngwriter -lpng `freetype-config --libs` -lfreetype -L.
OBJECTS = main.o wavelet.o CleanData.o ERSP.o convolution.o mechanical.o Plot.o 

morlet: $(OBJECTS)
	$(CXX) -o $@ $^ $(LDLIBS) $(CFLAGS);
	install_name_tool -change libpngwriter.dylib /opt/local/lib/libpngwriter.dylib morlet

clean:
	rm *.o morlet  *.log *.png