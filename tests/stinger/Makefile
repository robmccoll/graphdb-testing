include make.inc

#CORE
STINGER_CORE	= stinger.c stinger-deprecated.c stinger-iterator.c stinger-physmap.c stinger-return.c stinger-vertex.c stinger-workflow.c 
STINGER_CORE_SRC= $(addprefix src/core/, $(STINGER_CORE))
STINGER_CORE_OBJ= $(subst src,obj,$(subst .c,.o,$(STINGER_CORE_SRC)))

#UTIL
STINGER_UTIL	= csv.c histogram.c stinger-utils.c timer.c x86-full-empty.c xmalloc.c x86-full-empty.c graph_construct.c
STINGER_UTIL_SRC= $(addprefix src/util/, $(STINGER_UTIL))
STINGER_UTIL_OBJ= $(subst src,obj,$(subst .c,.o,$(STINGER_UTIL_SRC)))

#ALG
STINGER_ALG	= static_components.c static_multicontract_clustering.c streaming_components.c static_betweenness_centrality.c static_kcore.c streaming_clustering_coefficients.c result_writer.c static_pagerank.c
STINGER_ALG_SRC	= $(addprefix src/alg/, $(STINGER_ALG))
STINGER_ALG_OBJ	= $(subst src,obj,$(subst .c,.o,$(STINGER_ALG_SRC)))

#STREAM
STINGER_STREAM		= random_stream.c csv_stream.c binary_stream.c
STINGER_STREAM_SRC	= $(addprefix src/stream/, $(STINGER_STREAM))
STINGER_STREAM_OBJ	= $(subst src,obj,$(subst .c,.o,$(STINGER_STREAM_SRC)))

#LIB
STINGER_LIB	= mongoose/mongoose.c int-hm-seq/src/int-hm-seq.c string/src/astring.c fmemopen/src/fmemopen.c int-ht-seq/src/int-ht-seq.c
STINGER_LIB_SRC	= $(addprefix lib/,$(STINGER_LIB))
STINGER_LIB_OBJ	= $(subst src,obj,$(subst .c,.o,$(STINGER_LIB_SRC)))
STINGER_LIB_INCLUDE = -Ilib/mongoose -Ilib/string/inc -Ilib/fmemopen/inc -Ilib/int-ht-seq/inc -Ilib/int-hm-seq/inc

#FRAGMENTS - text files to be embedded as a string via a header
FRAGMENT     = vquery.html
FRAGMENT_SRC = $(addprefix include/fragments/,$(FRAGMENT))
FRAGMENT_OBJ = $(addsuffix .h,$(FRAGMENT_SRC))

#ALL
STINGER_ALL_SRC	= $(STINGER_CORE_SRC) $(STINGER_UTIL_SRC) $(STINGER_ALG_SRC) $(STINGER_STREAM_SRC) $(STINGER_LIB_SRC)
STINGER_ALL_OBJ	= $(STINGER_CORE_OBJ) $(STINGER_UTIL_OBJ) $(STINGER_ALG_OBJ) $(STINGER_STREAM_OBJ) $(STINGER_LIB_OBJ)

CFLAGS+= -Iinclude/alg -Iinclude/stream -Iinclude/util -Iinclude/core -Iinclude -I./ $(STINGER_LIB_INCLUDE)


include/fragments/%.h: include/fragments/%
	head -n 1 $^ > $@
	tail -n +2 $^ | sed -e 's/\\/\\\\/g' | sed -e 's/"/\\"/g' | sed -e 's/^/" /' | sed -e 's/$$/ \\n"/' >> $@
	echo ";" >> $@

.PHONY:	all
all: main 

printsettings:
	echo STINGER_CORE      $(STINGER_CORE)    
	echo STINGER_CORE_SRC  $(STINGER_CORE_SRC)
	echo STINGER_CORE_OBJ  $(STINGER_CORE_OBJ)
	echo STINGER_UTIL      $(STINGER_UTIL)    
	echo STINGER_UTIL_SRC  $(STINGER_UTIL_SRC)
	echo STINGER_UTIL_OBJ  $(STINGER_UTIL_OBJ)
	echo STINGER_ALG       $(STINGER_ALG)     
	echo STINGER_ALG_SRC   $(STINGER_ALG_SRC) 
	echo STINGER_ALG_OBJ   $(STINGER_ALG_OBJ) 
	echo STINGER_STREAM       $(STINGER_STREAM)     
	echo STINGER_STREAM_SRC   $(STINGER_STREAM_SRC) 
	echo STINGER_STREAM_OBJ   $(STINGER_STREAM_OBJ) 
	echo STINGER_ALL_SRC   $(STINGER_ALL_SRC) 
	echo STINGER_ALL_OBJ   $(STINGER_ALL_OBJ) 

objdirs:
	@ mkdir -p obj/core
	@ mkdir -p obj/util
	@ mkdir -p obj/alg
	@ mkdir -p obj/stream

lib/%:
	cd `echo $@ | sed -e 's/\(lib\/[^\/]*\)\/.*$$/\1/'`; make

main:	main.c $(STINGER_ALL_OBJ) $(BLECHIO)
	$(CC) $(MAINPLFLAG) $(CPPFLAGS) $(CFLAGS) -o $@ $^ \
		$(LDFLAGS) $(LDLIBS)

mainless:	main.c $(STINGER_ALL_OBJ) $(BLECHIO)
	$(CC) $(MAINPLFLAG) $(CPPFLAGS) $(CFLAGS) -o $@ $^ \
		$(LDFLAGS) $(LDLIBS) 2>&1 | less

lib: 	$(STINGER_ALL_OBJ)
	ar -cvq libstinger.a $^

obj/util/x86-full-empty.o:	src/util/x86-full-empty.c
	$(CC) $(MAINPLFLAG) $(CPPFLAGS) $(subst O1,O0,$(subst O2,O0,$(subst O3,O0,$(CFLAGS)))) -c -o $@ $^

static_example:	main.c r-utils.c src/timer.c src/xmalloc.c src/stinger-iterator.c $(BLECHIO) obj/x86-full-empty.o lib
	$(CC)  $(CFLAGS) -o static_main main.c libstinger.a 


lib/int-ht-seq/obj/%.o: lib/int-ht-seq/src/%.c
	$(CC) $(MAINPLFLAG) $(CPPFLAGS) $(CFLAGS) -fPIC -c $^ -o $@ 

obj/%.o: src/%.c | objdirs
	$(CC) $(MAINPLFLAG) $(CPPFLAGS) $(CFLAGS) -fPIC -c $^ -o $@ 


.PHONY: release
release:
	svn co http://svn.cc.gatech.edu/graphs/stinger/trunk trunk
	echo stinger-r`svn info trunk | grep "Revision: " | sed 's/.*Revision: //'` > RELEASE_DIR
	mv trunk `cat RELEASE_DIR`
	tar -czf `cat RELEASE_DIR`.tar.gz `cat RELEASE_DIR`
	rm -rf `cat RELEASE_DIR` RELEASE_DIR

.PHONY:	clean
clean:
	rm -rf main gen-streams workflow $(EXAMPLES) $(BLECHIO) $(BLECHIOGEN) \
		$(MAINPL) $(GENSTREAMSPL) libstinger.a $(STINGER_ALL_OBJ) include/fragments/*.h `find . -name "*.dSYM"`

######################################
# TESTS                              #
######################################
.PHONY: tests
tests: test.vertices test.physmap

test.vertices: obj/x86-full-empty.o src/stinger-vertex.c
	$(CC) $(MAINPLFLAG) $(CPPFLAGS) $(CFLAGS) -DSTINGER_VERTEX_TEST -o $@ $^ \
		$(LDFLAGS) $(LDLIBS)

test.physmap: src/stinger-physmap.c obj/x86-full-empty.o
	$(CC) $(MAINPLFLAG) -D PHYSMAP_TEST -Iinclude -o $@ $^
