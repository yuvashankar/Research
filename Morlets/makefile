CC      = gcc
CFLAGS  = -Wall -g -O2 -std=c99 -I../include

LDLIBS  = -lfftw3 -lgsl -lgslcblas -lm -fopenmp
OBJECTS = main.o wavelet.o CleanData.o ERSP.o

morlet: $(OBJECTS)
	$(CC) -o $@  $(OBJECTS) $(LDLIBS)

clean:
	rm *.o morlet  *.log
