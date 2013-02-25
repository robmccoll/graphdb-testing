.PHONY: all
all: main

main: test.c ../../lib/sqlite-amalgamation-3071502/sqlite3.c
	gcc -g -std=c99 -I ../../lib/sqlite-amalgamation-3071502/ -o $@ $^ -ldl -lpthread
