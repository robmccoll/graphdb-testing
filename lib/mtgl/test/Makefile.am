# -*- Makefile -*-

# traceview flags (see ../configure.in)
if WITH_TRACEVIEW
  TVFLAGS = -trace -trace_level 0 -profile
else
  TVFLAGS = 
endif

if WITH_DEBUG_ARRAYBOUNDS
  CCFLAG_ARRAYBOUNDS = -DCHECK_BOUNDS
else
  CCFLAG_ARRAYBOUNDS = 
endif

if ON_MTA
mostlyclean-local:
	/bin/rm -f $(EXTRA_PROGRAMS)
	/bin/rm -f *.ti
	/bin/rm -f *.ii
	/bin/rm -f *.pl
else
if ON_MAC
mostlyclean-local:
	/bin/rm -f $(EXTRA_PROGRAMS)
	/bin/rm -rf *.dSYM
else
mostlyclean-local:
	/bin/rm -f $(EXTRA_PROGRAMS)
endif
endif

# These programs are compiled via 'make check'.
check_PROGRAMS = \
	custom2file \
	graph2file \
	persist_graph \
	test_badrank \
	test_breadth_first_search \
	test_connected_components \
	test_expand_one_edge \
	test_independent_set \
	test_pagerank \
	test_pseudo_diameter \
	test_triangles \
	test_wcnm

# These programs are compiled individually via 'make <progname>'.
EXTRA_PROGRAMS = \
	generate_euler_tour \
	rmat2file \
	subiso_5cycles \
	subiso_triangles \
	switch_endianness \
	test_adjacency_iteration \
	test_duplicate \
	test_dynamic_array \
	test_edmonds_karp \
	test_edmonds_karp_small \
	test_filtering \
	test_graph \
	test_graph_parallel \
	test_hachar \
	test_mmap \
	test_mtgl_strtok \
	test_mutable_graph \
	test_numeric \
	test_partitioning_schemes \
	test_persist_xht \
	test_property_map \
	test_random \
	test_read_graph \
	test_shared_ptr \
	test_sort \
	test_sssp_deltastepping \
	test_subgraph \
	test_subgraph_isomorphism \
	test_vertex_betweenness \
	test_write_graph \
	test_xmt_hash_set \
	test_xmt_hash_table

buildtests: $(check_PROGRAMS)

buildall: $(check_PROGRAMS) $(EXTRA_PROGRAMS)

# Add these back in later when the test script is written.
#TESTS = $(noinst_PROGRAMS)
#XFAIL_TESTS = $(noinst_PROGRAMS)

# For check_PROGRAMS.
test_badrank_SOURCES = test_badrank.cpp
if ON_MTA
test_badrank_CFLAGS = -pl test_badrank.pl
test_badrank_CXXFLAGS = -pl test_badrank.pl
endif

test_breadth_first_search_SOURCES = test_breadth_first_search.cpp
if ON_MTA
test_breadth_first_search_CFLAGS = -pl test_breadth_first_search.pl
test_breadth_first_search_CXXFLAGS = -pl test_breadth_first_search.pl
endif

test_connected_components_SOURCES = test_connected_components.cpp
if ON_MTA
test_connected_components_CFLAGS = -pl test_connected_components.pl
test_connected_components_CXXFLAGS = -pl test_connected_components.pl
endif

test_expand_one_edge_SOURCES = test_expand_one_edge.cpp
if ON_MTA
test_expand_one_edge_CFLAGS = -pl test_expand_one_edge.pl
test_expand_one_edge_CXXFLAGS = -pl test_expand_one_edge.pl
endif

test_independent_set_SOURCES = test_independent_set.cpp
if ON_MTA
test_independent_set_CFLAGS = -pl test_independent_set.pl
test_independent_set_CXXFLAGS = -pl test_independent_set.pl
endif

test_pagerank_SOURCES = test_pagerank.cpp
if ON_MTA
test_pagerank_CFLAGS = -pl test_pagerank.pl
test_pagerank_CXXFLAGS = -pl test_pagerank.pl
endif

test_pseudo_diameter_SOURCES = test_pseudo_diameter.cpp
if ON_MTA
test_pseudo_diameter_CFLAGS = -pl test_pseudo_diameter.pl
test_pseudo_diameter_CXXFLAGS = -pl test_pseudo_diameter.pl
endif

test_triangles_SOURCES = test_triangles.cpp
if ON_MTA
test_triangles_CFLAGS = -pl test_triangles.pl
test_triangles_CXXFLAGS = -pl test_triangles.pl
endif

test_wcnm_SOURCES = test_wcnm.cpp ufl.cpp VolVolume.cpp
if ON_MTA
test_wcnm_CFLAGS = -pl test_wcnm.pl
test_wcnm_CXXFLAGS = -pl test_wcnm.pl
endif

custom2file_SOURCES = custom2file.cpp
if ON_MTA
custom2file_CFLAGS = -pl custom2file.pl
custom2file_CXXFLAGS = -pl custom2file.pl
endif

graph2file_SOURCES = graph2file.cpp
if ON_MTA
graph2file_CFLAGS = -pl graph2file.pl
graph2file_CXXFLAGS = -pl graph2file.pl
endif

persist_graph_SOURCES = persist_graph.cpp
if ON_MTA
persist_graph_CFLAGS = -pl persist_graph.pl
persist_graph_CXXFLAGS = -pl persist_graph.pl
endif


# For EXTRA_PROGRAMS.
generate_euler_tour_SOURCES = generate_euler_tour.cpp
if ON_MTA
generate_euler_tour_CFLAGS = -pl generate_euler_tour.pl
generate_euler_tour_CXXFLAGS = -pl generate_euler_tour.pl
endif

rmat2file_SOURCES = rmat2file.cpp
if ON_MTA
rmat2file_CFLAGS = -pl rmat2file.pl
rmat2file_CXXFLAGS = -pl rmat2file.pl
endif

subiso_5cycles_SOURCES = subiso_5cycles.cpp
if ON_MTA
subiso_5cycles_CFLAGS = -pl subiso_5cycles.pl
subiso_5cycles_CXXFLAGS = -pl subiso_5cycles.pl
endif

subiso_triangles_SOURCES = subiso_triangles.cpp
if ON_MTA
subiso_triangles_CFLAGS = -pl subiso_triangles.pl
subiso_triangles_CXXFLAGS = -pl subiso_triangles.pl
endif

switch_endianness_SOURCES = switch_endianness.cpp
if ON_MTA
switch_endianness_CFLAGS = -pl switch_endianness.pl
switch_endianness_CXXFLAGS = -pl switch_endianness.pl
endif

test_adjacency_iteration_SOURCES = test_adjacency_iteration.cpp
if ON_MTA
test_adjacency_iteration_CFLAGS = -pl test_adjacency_iteration.pl
test_adjacency_iteration_CXXFLAGS = -pl test_adjacency_iteration.pl
endif

test_duplicate_SOURCES = test_duplicate.cpp
if ON_MTA
test_duplicate_CFLAGS = -pl test_duplicate.pl
test_duplicate_CXXFLAGS = -pl test_duplicate.pl
endif

test_dynamic_array_SOURCES = test_dynamic_array.cpp
if ON_MTA
test_dynamic_array_CFLAGS = -pl test_dynamic_array.pl
test_dynamic_array_CXXFLAGS = -pl test_dynamic_array.pl
endif

test_edmonds_karp_SOURCES = test_edmonds_karp.cpp
if ON_MTA
test_edmonds_karp_CFLAGS = -pl test_edmonds_karp.pl
test_edmonds_karp_CXXFLAGS = -pl test_edmonds_karp.pl
endif

test_edmonds_karp_small_SOURCES = test_edmonds_karp_small.cpp
if ON_MTA
test_edmonds_karp_small_CFLAGS = -pl test_edmonds_karp_small.pl
test_edmonds_karp_small_CXXFLAGS = -pl test_edmonds_karp_small.pl
endif

test_filtering_SOURCES = test_filtering.cpp
if ON_MTA
test_filtering_CFLAGS = -pl test_filtering.pl
test_filtering_CXXFLAGS = -pl test_filtering.pl
endif

test_graph_SOURCES = test_graph.cpp
if ON_MTA
test_graph_CFLAGS = -pl test_graph.pl
test_graph_CXXFLAGS = -pl test_graph.pl
endif

test_graph_parallel_SOURCES = test_graph_parallel.cpp
if ON_MTA
test_graph_parallel_CFLAGS = -pl test_graph_parallel.pl
test_graph_parallel_CXXFLAGS = -pl test_graph_parallel.pl
endif

test_hachar_SOURCES = test_hachar.cpp
if ON_MTA
test_hachar_CFLAGS = -pl test_hachar.pl
test_hachar_CXXFLAGS = -pl test_hachar.pl
endif

test_mmap_SOURCES = test_mmap.cpp
if ON_MTA
test_mmap_CFLAGS = -pl test_mmap.pl
test_mmap_CXXFLAGS = -pl test_mmap.pl
endif

test_mtgl_strtok_SOURCES = test_mtgl_strtok.cpp
if ON_MTA
test_mtgl_strtok_CFLAGS = -pl test_mtgl_strtok.pl
test_mtgl_strtok_CXXFLAGS = -pl test_mtgl_strtok.pl
endif

test_mutable_graph_SOURCES = test_mutable_graph.cpp
if ON_MTA
test_mutable_graph_CFLAGS = -pl test_mutable_graph.pl
test_mutable_graph_CXXFLAGS = -pl test_mutable_graph.pl
endif

test_numeric_SOURCES = test_numeric.cpp
if ON_MTA
test_numeric_CFLAGS = -pl test_numeric.pl
test_numeric_CXXFLAGS = -pl test_numeric.pl
endif

test_partitioning_schemes_SOURCES = test_partitioning_schemes.cpp
if ON_MTA
test_partitioning_schemes_CFLAGS = -pl test_partitioning_schemes.pl
test_partitioning_schemes_CXXFLAGS = -pl test_partitioning_schemes.pl
endif

test_persist_xht_SOURCES = test_persist_xht.cpp
if ON_MTA
test_persist_xht_CFLAGS = -pl test_persist_xht.pl
test_persist_xht_CXXFLAGS = -pl test_persist_xht.pl
endif

test_property_map_SOURCES = test_property_map.cpp
if ON_MTA
test_property_map_CFLAGS = -pl test_property_map.pl
test_property_map_CXXFLAGS = -pl test_property_map.pl
endif

test_random_SOURCES = test_random.cpp
if ON_MTA
test_random_CFLAGS = -pl test_random.pl
test_random_CXXFLAGS = -pl test_random.pl
endif

test_read_graph_SOURCES = test_read_graph.cpp
if ON_MTA
test_read_graph_CFLAGS = -pl test_read_graph.pl
test_read_graph_CXXFLAGS = -pl test_read_graph.pl
endif

test_shared_ptr_SOURCES = test_shared_ptr.cpp
if ON_MTA
test_shared_ptr_CFLAGS = -pl test_shared_ptr.pl
test_shared_ptr_CXXFLAGS = -pl test_shared_ptr.pl
endif

test_sort_SOURCES = test_sort.cpp
if ON_MTA
test_sort_CFLAGS = -pl test_sort.pl
test_sort_CXXFLAGS = -pl test_sort.pl
endif

test_sssp_deltastepping_SOURCES = test_sssp_deltastepping.cpp
if ON_MTA
test_sssp_deltastepping_CFLAGS = -pl test_sssp_deltastepping.pl
test_sssp_deltastepping_CXXFLAGS = -pl test_sssp_deltastepping.pl
endif

test_subgraph_SOURCES = test_subgraph.cpp
if ON_MTA
test_subgraph_CFLAGS = -pl test_subgraph.pl
test_subgraph_CXXFLAGS = -pl test_subgraph.pl
endif

test_subgraph_isomorphism_SOURCES = test_subgraph_isomorphism.cpp
if ON_MTA
test_subgraph_isomorphism_CFLAGS = -pl test_subgraph_isomorphism.pl
test_subgraph_isomorphism_CXXFLAGS = -pl test_subgraph_isomorphism.pl
endif

test_vertex_betweenness_SOURCES = test_vertex_betweenness.cpp
if ON_MTA
test_vertex_betweenness_CFLAGS = -pl test_vertex_betweenness.pl
test_vertex_betweenness_CXXFLAGS = -pl test_vertex_betweenness.pl
endif

test_write_graph_SOURCES = test_write_graph.cpp
if ON_MTA
test_write_graph_CFLAGS = -pl test_write_graph.pl
test_write_graph_CXXFLAGS = -pl test_write_graph.pl
endif

test_xmt_hash_set_SOURCES = test_xmt_hash_set.cpp
if ON_MTA
test_xmt_hash_set_CFLAGS = -pl test_xmt_hash_set.pl
test_xmt_hash_set_CXXFLAGS = -pl test_xmt_hash_set.pl
endif

test_xmt_hash_table_SOURCES = test_xmt_hash_table.cpp
if ON_MTA
test_xmt_hash_table_CFLAGS = -pl test_xmt_hash_table.pl
test_xmt_hash_table_CXXFLAGS = -pl test_xmt_hash_table.pl
endif
