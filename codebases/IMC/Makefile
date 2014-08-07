CXX ?= g++
CC ?= gcc
CFLAGS = -w -fopenmp -Wconversion -O3 -fPIC
LIBS = blas/blas.a
SHVER = 1
OS = $(shell uname)
#LIBS = -lblas

all: clean predict train

train: tron.o  blas/blas.a util.o
	$(CXX) $(CFLAGS) -o train train.cpp tron.o util.o $(LIBS)

predict: util.o
	$(CXX) $(CFLAGS) -o predict predict.cpp util.o

tron.o: tron.cpp tron.h
	$(CXX) $(CFLAGS) -c -o tron.o tron.cpp

blas/blas.a: blas/*.c blas/*.h
	make -C blas OPTFLAGS='$(CFLAGS)' CC='$(CC)';
util.o:
	$(CXX) $(CFLAGS) -c -o util.o util.cpp
clean:
	make -C blas clean
	rm -f *~ *.o train predict
