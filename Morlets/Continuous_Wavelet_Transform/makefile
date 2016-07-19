CC      = gcc
CFLAGS  = -Wall -lfftw3 -lgsl -lgslcblas -lm -g -O2 -fopenmp
OBJECTS = main.o wavelet.o

cwt: $(OBJECTS)
	$(CC) -o $@  $(OBJECTS) $(CFLAGS)

clean:
	rm *.o cwt  *.log
	rm -rf *.dSYM
