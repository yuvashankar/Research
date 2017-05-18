CC  = gcc
IDIR =../include
CFLAGS= -I$(IDIR) -I/opt/local/lib/gcc5/gcc/x86_64-apple-darwin15/5.4.0/include -Wall -Wextra -g -fopenmp -pedantic-errors

LDLIBS  = -lfftw3 -lfftw3_threads -lgsl -lgslcblas -lm -L/opt/local/lib/

ODIR=obj
LDIR =../lib

_DEPS = edflib.h wavelet.h processEEG.h ERSP.h pngwriter.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o edflib.o wavelet.o ERSP.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

ERSP: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LDLIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ *.log

