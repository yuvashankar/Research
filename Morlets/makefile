CC      = gcc
CFLAGS  = -Wall -lfftw3 -lgsl -lgslcblas -lm -g -O2 -fopenmp -DNDEBUG
OBJECTS = main.o Morlet.o wavelet.o

morlet: $(OBJECTS)
	$(CC) -o $@  $(OBJECTS) $(CFLAGS)

clean:
	rm *.o morlet  *.log
