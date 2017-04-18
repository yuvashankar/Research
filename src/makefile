# CC  = gcc
# CXX = g++
# IDIR =../include
# CFLAGS= -I$(IDIR) -I/opt/local/lib/gcc5/gcc/x86_64-apple-darwin15/5.4.0/include -Wall -Wextra -g -fopenmp -pedantic-errors
# CXXFLAGS = -Wall -Wextra -g -fopenmp `freetype-config --cflags` -I/opt/local/lib/gcc5/gcc/x86_64-apple-darwin15/5.4.0/include -I../include

# LDLIBS  = -lfftw3 -lfftw3_threads -lgsl -lgslcblas -lm -lz `freetype-config --libs` -lz -lpngwriter -lpng

# ODIR=obj
# LDIR =../lib

# _DEPS = edflib.h wavelet.h processEEG.h ERSP.h pngwriter.h
# DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

# _OBJ = main.o edflib.o OpenFile.o FindTriggers.o FilterTriggers.o wavelet.o CleanData.o ERSP.o
# OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

# $(ODIR)/%.o: %.c $(DEPS)
# 	$(CXX) -c -o $@ $< $(CFLAGS)

# processEEG: $(OBJ)
# 	gcc -o $@ $^ $(CFLAGS) $(LDLIBS)
# 	install_name_tool -change libpngwriter.dylib /opt/local/lib/libpngwriter.dylib morlet

# .PHONY: clean

# clean:
# 	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ *.log

CC       = gcc
CXX      = g++
CFLAGS   = -I/opt/local/lib/gcc5/gcc/x86_64-apple-darwin15/5.4.0/include -I../include -Wall -Wextra -g -fopenmp -pedantic-errors
CXXFLAGS = -Wall -Wextra -g -fopenmp `freetype-config --cflags` -I/opt/local/lib/gcc5/gcc/x86_64-apple-darwin15/5.4.0/include -I../include

LDLIBS  = -lfftw3 -lfftw3_threads -lgsl -lgslcblas -lm -lz `freetype-config --libs` -lz -lpngwriter -lpng
OBJECTS = main.o edflib.o OpenFile.o FindTriggers.o FilterTriggers.o wavelet.o CleanData.o ERSP.o

morlet: $(OBJECTS)
	$(CXX) -o $@ $^ $(CFLAGS) $(LDLIBS);
	install_name_tool -change libpngwriter.dylib /opt/local/lib/libpngwriter.dylib morlet

clean:
	rm *.o morlet  *.log *.png
