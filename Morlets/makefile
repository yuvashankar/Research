CC      = gcc
CFLAGS  = -Wall -lfftw3 -lgsl -lgslcblas -lm
OBJECTS = main.o Morlet.o wavelet.o

morlet: $(OBJECTS)
	$(CC) -o $@  $(OBJECTS) $(CFLAGS)

clean:
	rm *.o morlet  *.log
