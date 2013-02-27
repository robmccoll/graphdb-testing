LIB_DIR=../../lib/
LIB_INCLUDES=timer dexcpp-4.7.1/includes/stlport dexcpp-4.7.1/includes/dex
LIB=-L$(LIB_DIR)/dexcpp-4.7.1/lib/linux64 -ldex -lstlport -lpthread

INCLUDE=$(addprefix -I$(LIB_DIR), $(LIB_INCLUDES))

DEFINE=-D_FILE_OFFSET_BITS=64


.PHONY: all
all: main

timer.o: ../../lib/timer/timer.c
	gcc -g -c -O2 $(INCLUDE) $(DEFINE) -o $@ $^ -lm -lrt

main: test.cpp timer.o
	g++ -g -O2 -std=c++0x $(INCLUDE) $(DEFINE) -o $@ $^ -lm -lrt $(LIB)
