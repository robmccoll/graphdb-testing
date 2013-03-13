CC=gcc -std=gnu9x
CXX=g++
CFLAGS=-fopenmp -g -O2
LDLIBS=-lm -lrt

CFLAGS+=-Iinc

rmatter: src/rmatter.c src/random.c src/timer.c
	$(CC) $(MAINPLFLAG) $(CPPFLAGS) $(CFLAGS) -o $@ $^ \
		$(LDFLAGS) $(LDLIBS)

.PHONY: all
all: rmatter

.PHONY:	clean
clean:
	rm -f main gen-streams rmatter $(EXAMPLES) $(BLECHIO) $(BLECHIOGEN) \
		$(MAINPL) $(GENSTREAMSPL) libstinger.a obj/*.o

