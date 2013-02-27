# -*- Makefile -*-

# These are extra files that need to be distributed.
EXTRA_DIST = \
	graphs \
	subgraph_isomorphism1.cpp \
	subgraph_isomorphism2.cpp \
	subgraph_isomorphism3.cpp \
	subgraph_isomorphism4.cpp

if ON_MTA
mostlyclean-local:
	/bin/rm -f $(EXTRA_PROGRAMS)
	/bin/rm -f *.ti
	/bin/rm -f *.ii
	/bin/rm -f *.pl
else
mostlyclean-local:
	/bin/rm -f $(EXTRA_PROGRAMS)
endif

# These programs are compiled via 'make tutorial'.
EXTRA_PROGRAMS = \
	example1 \
	example2 \
	example3 \
	exercise1 \
	exercise2 \
	exercise3 \
	exercise4 \
	filter_example \
	input_output \
	property_maps1 \
	property_maps2 \
	property_maps3 \
	property_maps4 \
	property_maps5 \
	triangles

tutorial: $(EXTRA_PROGRAMS)

.PHONY: tutorial

# For EXTRA_PROGRAMS.
example1_SOURCES = example1.cpp
if ON_MTA
example1_CFLAGS = -pl example1.pl
example1_CXXFLAGS = -pl example1.pl
endif

example2_SOURCES = example2.cpp
if ON_MTA
example2_CFLAGS = -pl example2.pl
example2_CXXFLAGS = -pl example2.pl
endif

example3_SOURCES = example3.cpp
if ON_MTA
example3_CFLAGS = -pl example3.pl
example3_CXXFLAGS = -pl example3.pl
endif

exercise1_SOURCES = exercise1.cpp
if ON_MTA
exercise1_CFLAGS = -pl exercise1.pl
exercise1_CXXFLAGS = -pl exercise1.pl
endif

exercise2_SOURCES = exercise2.cpp
if ON_MTA
exercise2_CFLAGS = -pl exercise2.pl
exercise2_CXXFLAGS = -pl exercise2.pl
endif

exercise3_SOURCES = exercise3.cpp
if ON_MTA
exercise3_CFLAGS = -pl exercise3.pl
exercise3_CXXFLAGS = -pl exercise3.pl
endif

exercise4_SOURCES = exercise4.cpp
if ON_MTA
exercise4_CFLAGS = -pl exercise4.pl
exercise4_CXXFLAGS = -pl exercise4.pl
endif

filter_example_SOURCES = filter_example.cpp
if ON_MTA
filter_example_CFLAGS = -pl filter_example.pl
filter_example_CXXFLAGS = -pl filter_example.pl
endif

input_output_SOURCES = input_output.cpp
if ON_MTA
input_output_CFLAGS = -pl input_output.pl
input_output_CXXFLAGS = -pl input_output.pl
endif

property_maps1_SOURCES = property_maps1.cpp
if ON_MTA
property_maps1_CFLAGS = -pl property_maps1.pl
property_maps1_CXXFLAGS = -pl property_maps1.pl
endif

property_maps2_SOURCES = property_maps2.cpp
if ON_MTA
property_maps2_CFLAGS = -pl property_maps2.pl
property_maps2_CXXFLAGS = -pl property_maps2.pl
endif

property_maps3_SOURCES = property_maps3.cpp
if ON_MTA
property_maps3_CFLAGS = -pl property_maps3.pl
property_maps3_CXXFLAGS = -pl property_maps3.pl
endif

property_maps4_SOURCES = property_maps4.cpp
if ON_MTA
property_maps4_CFLAGS = -pl property_maps4.pl
property_maps4_CXXFLAGS = -pl property_maps4.pl
endif

property_maps5_SOURCES = property_maps5.cpp
if ON_MTA
property_maps5_CFLAGS = -pl property_maps5.pl
property_maps5_CXXFLAGS = -pl property_maps5.pl
endif

triangles_SOURCES = triangles.cpp
if ON_MTA
triangles_CFLAGS = -pl triangles.pl
triangles_CXXFLAGS = -pl triangles.pl
endif
