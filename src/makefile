IDIR =../include
CC=gcc
CFLAGS=-I$(IDIR) -Wall -lfftw3 -lgsl -fopenmp -g -O2

ODIR=obj
LDIR =../lib

_DEPS = edflib.h wavelet.h processEEG.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o edflib.o OpenFile.o FindTriggers.o FilterTriggers.o wavelet.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

processEEG: $(OBJ)
	gcc -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ *.log
