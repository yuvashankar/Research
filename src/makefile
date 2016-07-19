IDIR =../include
CC=gcc
CFLAGS=-I$(IDIR) -Wall -g -O2 -std=c99 -I../../bin/include -I/usr/include/openmpi-x86_64 -I/usr/local/cuda/include

LDFLAGS = -L../../bin/lib
LDLIBS = -lfftw3 -lgsl -lgslcblas -lm -fopenmp

ODIR=obj
LDIR =../lib

_DEPS = edflib.h wavelet.h processEEG.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o edflib.o OpenFile.o FindTriggers.o FilterTriggers.o wavelet.o CleanData.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

processEEG: $(OBJ)
	gcc -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ *.log
