.PHONY: all
all: main

timer.o: ../../lib/timer/timer.c
	gcc -I ../../lib/timer -c -o $@ $^ 

main: test.cpp timer.o
	g++ -g -O2 -I ../../lib/timer -I ../../lib/boost -o $@ $^ -lm -lrt
