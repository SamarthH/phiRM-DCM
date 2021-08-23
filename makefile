CC=gcc
CFLAGS= -lm -lgsl -lgslcblas -Wall -fopenmp -Wl,-R/usr/local/lib
DEPS = phase_fit.h
OBJ = phiRM_DCM.o phase_fit.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

resonanceFactor: $(OBJ)
	$(CC) -o resonanceFactor $(OBJ) $(CFLAGS)

clean:
	rm *.o
