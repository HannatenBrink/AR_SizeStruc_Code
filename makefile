C=g++
CFLAGS= -std=c++11 -march=native -Wall -O3

all: RunIBM

RunIBM: main.o Individual.o 
	$(C) -L/usr/local/include main.o  Individual.o   -o RunIBM

main.o: main.cpp main.h Individual.h
	$(C) $(CFLAGS) -c main.cpp -o main.o


Individual.o: Individual.cpp Individual.h
	$(C) $(CFLAGS) -c Individual.cpp -o Individual.o


clean:
	rm *o RunIBM
