IDIR =../include
CC=gcc
CFLAGS=-I$(IDIR) -Wall -g -O2 -std=c99 -fopenmp

LDLIBS = -lfftw3 -lgsl -lgslcblas -lm

ODIR=obj
LDIR =../lib

_DEPS = edflib.h wavelet.h processEEG.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o edflib.o OpenFile.o FindTriggers.o FilterTriggers.o wavelet.o CleanData.o ERSP.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

processEEG: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LDLIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ *.log
