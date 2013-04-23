#CONFIG
APPNAME=int-hm-seq
CSRC=int-hm-seq.c timer.c
IS_LIB=yes

OBJDIR=obj
SRCDIR=src

CC=gcc
INCLUDES=./ ./inc
CFLAGS=-std=gnu99 -O3 -g
LIBS=

include ../../tools/default.make

$(APPNAME)-test: $(addprefix $(SRCDIR)/,$(CSRC))
	$(CC) $(CFLAGS) -DINT_HM_SEQ_TEST $^ -o $@ $(LIBS) -lrt
