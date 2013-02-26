.PHONY: all
all: main

timer.o: ../../lib/timer/timer.c
	gcc -I ../../lib/timer -c -o $@ $^ 

main: test.cpp timer.o
	g++ -g -fopenmp -O2 -I ../../lib/timer -I ../../lib/mtgl/build/include -o $@ $^ -lm -lrt
