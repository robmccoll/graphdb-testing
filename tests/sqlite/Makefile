.PHONY: all
all: main

timer.o: ../../lib/timer/timer.c
	gcc -I ../../lib/timer -c -o $@ $^ 

main: timer.o test.c ../../lib/sqlite-amalgamation-3071502/sqlite3.c
	gcc -g -std=c99 -I ../../lib/timer -I ../../lib/sqlite-amalgamation-3071502/ -o $@ $^ -ldl -lpthread -lrt -lm
