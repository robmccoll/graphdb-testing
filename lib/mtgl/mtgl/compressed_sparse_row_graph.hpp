/*  _________________________________________________________________________
 *
 *  MTGL: The MultiThreaded Graph Library
 *  Copyright (c) 2008 Sandia Corporation.
 *  This software is distributed under the BSD License.
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *  For more information, see the README file in the top MTGL directory.
 *  _________________________________________________________________________
 */

/****************************************************************************/
/*! \file compressed_sparse_row_graph.hpp

    \brief This graph class implements an enhanced compressed sparse row
           graph that can be directed, undirected, or bidirectional.

    \author Greg Mackey (gemacke@sandia.gov)
    \author Jon Berry (jberry@sandia.gov)
    \author Brad Mancke
    \author Kamesh Madduri

    \date 1/3/2008

    The first enhancement is storing the source points of the adjacencies in
    addition to the end points.  This yields direct access to source vertices
    of an edge without the need for searching through the index array (aka it
    runs faster). The second enhancement is assigning edge ids based on the
    order the edges are passed during the call to init() instead of based on
    their indices.  These edge ids are kept as two arrays that allow direct
    access both directions between the edge id and the edge index.  This
    allows users to create array property maps for edges more easily as a user
    usually has the property array in the same order as they submitted the
    edges to init().
*/
/****************************************************************************/

#ifndef MTGL_COMPRESSED_SPARSE_ROW_GRAPH_HPP
#define MTGL_COMPRESSED_SPARSE_ROW_GRAPH_HPP

#include <limits>

#include <mtgl/mtgl_adapter.hpp>

#ifdef USING_QTHREADS
#include <qthread/qloop.hpp>
#endif

namespace mtgl {

namespace detail {

#ifdef USING_QTHREADS
template <typename T>
class increment_degrees {
public:
  increment_degrees(T* d, T* v) : degree(d), verts(v) {}

  void operator()(const size_t start, const size_t stop)
  {
    for (size_t i = start; i != stop; ++i) mt_incr(degree[verts[i]], 1);
  }

private:
  T* degree;
  T* verts;
};

template <typename T>
class reset_degrees {
public:
  reset_degrees(T* d) : degree(d) {}

  void operator()(const size_t start, const size_t stop)
  {
    for (size_t i = start; i != stop; ++i) degree[i] = 0;
  }

private:
  T* degree;
};

template <typename T>
class add_csr_edges {
public:
  add_csr_edges(T* dg, T* s, T* d, T* ne, T* sv, T* ev, T* lo, T* li) :
    degree(dg), srcs(s), dests(d), numEdges(ne), srcV(sv), endV(ev),
    loc_original_ids(lo), loc_internal_ids(li) {}

  void operator()(const size_t start, const size_t stop)
  {
    for (size_t i = start; i != stop; ++i)
    {
      T u = srcs[i];
      T vpos = numEdges[u] + mt_incr(degree[u], 1);
      srcV[vpos] = srcs[i];
      endV[vpos] = dests[i];
      loc_original_ids[vpos] = i;
      loc_internal_ids[i] = vpos;
    }
  }

private:
  T* degree;
  T* srcs;
  T* dests;
  T* numEdges;
  T* srcV;
  T* endV;
  T* loc_original_ids;
  T* loc_internal_ids;
};

template <typename T>
class add_csr_rev_edges {
public:
  add_csr_rev_edges(T* dg, T* s, T* d, T* ne, T* sv, T* ev, T* lo) :
    degree(dg), srcs(s), dests(d), numEdges(ne), srcV(sv), endV(ev),
    loc_original_ids(lo) {}

  void operator()(const size_t start, const size_t stop)
  {
    for (size_t i = start; i != stop; ++i)
    {
      T u = dests[i];
      T vpos = numEdges[u] + mt_incr(degree[u], 1);
      srcV[vpos] = dests[i];
      endV[vpos] = srcs[i];
      loc_original_ids[vpos] = i;
    }
  }

private:
  T* degree;
  T* srcs;
  T* dests;
  T* numEdges;
  T* srcV;
  T* endV;
  T* loc_original_ids;
};

template <typename T>
class add_csr_bid_edges {
public:
  add_csr_bid_edges(T* od, T* id, T* s, T* d, T* one, T* ine, T* osv, T* oev,
                    T* isv, T* iev, T* lo, T* li, T* cci) :
    out_degree(od), in_degree(id), srcs(s), dests(d), out_numEdges(one),
    in_numEdges(ine), out_srcV(osv), out_endV(oev), in_srcV(isv), in_endV(iev),
    loc_original_ids(lo), loc_internal_ids(li), cross_cross_index(cci) {}

  void operator()(const size_t start, const size_t stop)
  {
    for (size_t i = start; i != stop; ++i)
    {
      T src = srcs[i];
      T dest = dests[i];

      // Add the edge to the source vertex.
      T out_vp = out_numEdges[src] + mt_incr(out_degree[src], 1);
      out_endV[out_vp] = dest;
      out_srcV[out_vp] = src;
      loc_original_ids[out_vp] = i;
      loc_internal_ids[i] = out_vp;

      // Add the edge to the dest vertex.
      T in_vp = in_numEdges[dest] + mt_incr(in_degree[dest], 1);
      in_endV[in_vp] = src;
      in_srcV[in_vp] = dest;
      cross_cross_index[in_vp] = out_vp;
    }
  }

private:
  T* out_degree;
  T* in_degree;
  T* srcs;
  T* dests;
  T* out_numEdges;
  T* in_numEdges;
  T* out_srcV;
  T* out_endV;
  T* in_srcV;
  T* in_endV;
  T* loc_original_ids;
  T* loc_internal_ids;
  T* cross_cross_index;
};

template <typename T>
class add_rev_original_ids {
public:
  add_rev_original_ids(T* lo, T* lro, T* cci) :
    loc_original_ids(lo), loc_rev_original_ids(lro), cross_cross_index(cci) {}

  void operator()(const size_t start, const size_t stop)
  {
    for (size_t i = start; i != stop; ++i)
    {
      loc_rev_original_ids[i] = loc_original_ids[cross_cross_index[i]];
    }
  }

private:
  T* loc_original_ids;
  T* loc_rev_original_ids;
  T* cross_cross_index;
};

#ifndef MTGL_ARR_COPY_LOOP
#define MTGL_ARR_COPY_LOOP
template <typename T>
class arr_copy_loop {
public:
  arr_copy_loop(T* d, T* s) : dest(d), src(s) {}

  void operator()(const size_t start, const size_t stop)
  {
    for (size_t i = start; i != stop; ++i) dest[i] = src[i];
  }

private:
  T* dest;
  T* src;
};
#endif
#endif

template <typename Graph>
class csr_edge_adapter {
private:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;

public:
  csr_edge_adapter() : first((std::numeric_limits<vertex_descriptor>::max)()),
                       second((std::numeric_limits<vertex_descriptor>::max)()),
                       id((std::numeric_limits<size_type>::max)()) {}

  csr_edge_adapter(vertex_descriptor v1, vertex_descriptor v2, size_type eid) :
    first(v1), second(v2), id(eid) {}

  bool operator==(const csr_edge_adapter& rhs) const
  { return id == rhs.id; }

public:
  vertex_descriptor first;
  vertex_descriptor second;
  size_type id;
};

/***/

template <typename Graph>
class csr_vertex_iterator {
private:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;

public:
  csr_vertex_iterator() {}
  vertex_descriptor operator[](size_type p) const { return p; }
};

/***/

template <typename Graph>
class csr_edge_iterator {
private:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

public:
  csr_edge_iterator() : src_points(0), end_points(0), internal_ids(0) {}

  csr_edge_iterator(size_type* sp, size_type* ep, size_type* ii) :
    src_points(sp), end_points(ep), internal_ids(ii) {}

  edge_descriptor operator[](size_type p) const
  {
    return edge_descriptor(src_points[internal_ids[p]],
                           end_points[internal_ids[p]], p);
  }

private:
  size_type* src_points;
  size_type* end_points;
  size_type* internal_ids;
};

/***/

template <typename Graph>
class csr_adjacency_iterator {
private:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;

public:
  csr_adjacency_iterator() : adj(0) {}

  csr_adjacency_iterator(size_type* ind, size_type* adj_, size_type i) :
    adj(adj_ + ind[i]) {}

  vertex_descriptor operator[](size_type p) const { return adj[p]; }

private:
  size_type* adj;
};

/***/

template <typename Graph>
class csr_out_edge_iterator {
private:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

public:
  csr_out_edge_iterator() : adj(0), original_ids(0), vid(0) {}

  csr_out_edge_iterator(size_type* ind, size_type* adj_, size_type* oi,
                        size_type i) :
    adj(adj_ + ind[i]), original_ids(oi + ind[i]), vid(i) {}

  edge_descriptor operator[](size_type p) const
  {
    return edge_descriptor(vid, adj[p], original_ids[p]);
  }

private:
  size_type* adj;
  size_type* original_ids;
  size_type vid;
};

/***/

template <typename Graph>
class csr_thread_vertex_iterator {
private:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;

public:
  csr_thread_vertex_iterator() : pos((std::numeric_limits<size_type>::max)()) {}

  csr_thread_vertex_iterator(size_type v) : pos(v) {}

  csr_thread_vertex_iterator(const csr_thread_vertex_iterator& rhs) :
    pos(rhs.pos) {}

  csr_thread_vertex_iterator& operator=(const csr_thread_vertex_iterator& rhs)
  {
    pos = rhs.pos;
    return *this;
  }

  csr_thread_vertex_iterator& operator++()
  {
    ++pos;
    return *this;
  }

  csr_thread_vertex_iterator& operator++(int)
  {
    csr_thread_vertex_iterator temp(*this);

    ++pos;
    return temp;
  }

  vertex_descriptor operator*() const { return pos; }

  bool operator==(const csr_thread_vertex_iterator& rhs) const
  { return pos == rhs.pos; }
  bool operator!=(const csr_thread_vertex_iterator& rhs) const
  { return pos != rhs.pos; }

private:
  size_type pos;
};

/***/

template <typename Graph>
class csr_thread_edge_iterator {
private:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

public:
  csr_thread_edge_iterator() : pos((std::numeric_limits<size_type>::max)()),
                               src_points(0), end_points(0), internal_ids(0) {}

  csr_thread_edge_iterator(vertex_descriptor* sp, vertex_descriptor* ep,
                           size_type* ii, size_type e) :
    pos(e), src_points(sp), end_points(ep), internal_ids(ii) {}

  csr_thread_edge_iterator(const csr_thread_edge_iterator& rhs) :
    pos(rhs.pos), src_points(rhs.src_points), end_points(rhs.end_points),
    internal_ids(rhs.internal_ids) {}

  csr_thread_edge_iterator& operator=(const csr_thread_edge_iterator& rhs)
  {
    pos = rhs.pos;
    src_points = rhs.src_points;
    end_points = rhs.end_points;
    internal_ids = rhs.internal_ids;
    return *this;
  }

  csr_thread_edge_iterator& operator++()
  {
    ++pos;
    return *this;
  }

  csr_thread_edge_iterator& operator++(int)
  {
    csr_thread_edge_iterator temp(*this);

    ++pos;
    return temp;
  }

  edge_descriptor operator*() const
  {
    return edge_descriptor(src_points[internal_ids[pos]],
                           end_points[internal_ids[pos]], pos);
  }

  bool operator==(const csr_thread_edge_iterator& rhs) const
  { return pos == rhs.pos; }
  bool operator!=(const csr_thread_edge_iterator& rhs) const
  { return pos != rhs.pos; }

private:
  size_type pos;
  vertex_descriptor* src_points;
  vertex_descriptor* end_points;
  size_type* internal_ids;
};

/***/

template <typename Graph>
class csr_thread_adjacency_iterator {
private:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;

public:
  csr_thread_adjacency_iterator() : adj(0) {}

  csr_thread_adjacency_iterator(size_type* ind, vertex_descriptor* adj_,
                                size_type i, size_type adjp) :
    adj(adj_ + ind[i] + adjp) {}

  csr_thread_adjacency_iterator(const csr_thread_adjacency_iterator& rhs) :
    adj(rhs.adj) {}

  csr_thread_adjacency_iterator&
  operator=(const csr_thread_adjacency_iterator& rhs)
  {
    adj = rhs.adj;
    return *this;
  }

  csr_thread_adjacency_iterator& operator++()
  {
    ++adj;
    return *this;
  }

  csr_thread_adjacency_iterator& operator++(int)
  {
    csr_thread_adjacency_iterator temp(*this);

    ++adj;
    return temp;
  }

  vertex_descriptor operator*() const { return *adj; }

  bool operator==(const csr_thread_adjacency_iterator& rhs)
  { return adj == rhs.adj; }
  bool operator!=(const csr_thread_adjacency_iterator& rhs)
  { return adj != rhs.adj; }

private:
  vertex_descriptor* adj;
};

/***/

template <typename Graph>
class csr_thread_out_edge_iterator {
private:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

public:
  csr_thread_out_edge_iterator() : adj(0), original_ids(0), vid(0) {}

  csr_thread_out_edge_iterator(size_type* ind, vertex_descriptor* adj_,
                               size_type* oi, size_type i, size_type adjp) :
    adj(adj_ + ind[i] + adjp), original_ids(oi + ind[i] + adjp), vid(i) {}

  csr_thread_out_edge_iterator(const csr_thread_out_edge_iterator& rhs) :
    adj(rhs.adj), original_ids(rhs.original_ids), vid(rhs.vid) {}

  csr_thread_out_edge_iterator&
  operator=(const csr_thread_out_edge_iterator& rhs)
  {
    adj = rhs.adj;
    original_ids = rhs.original_ids;
    vid = rhs.vid;
    return *this;
  }

  csr_thread_out_edge_iterator& operator++()
  {
    ++adj;
    ++original_ids;
    return *this;
  }

  csr_thread_out_edge_iterator& operator++(int)
  {
    csr_thread_out_edge_iterator temp(*this);

    ++adj;
    ++original_ids;
    return temp;
  }

  edge_descriptor operator*() const
  {
    return edge_descriptor(vid, *adj, *original_ids);
  }

  bool operator==(const csr_thread_out_edge_iterator& rhs)
  { return adj == rhs.adj; }
  bool operator!=(const csr_thread_out_edge_iterator& rhs)
  { return adj != rhs.adj; }

private:
  vertex_descriptor* adj;
  size_type* original_ids;
  size_type vid;
};

}

/***/

template <typename DIRECTION = directedS>
class compressed_sparse_row_graph {
public:
  typedef unsigned long size_type;
  typedef unsigned long vertex_descriptor;
  typedef detail::csr_edge_adapter<compressed_sparse_row_graph> edge_descriptor;
  typedef detail::csr_vertex_iterator<compressed_sparse_row_graph>
          vertex_iterator;
  typedef detail::csr_adjacency_iterator<compressed_sparse_row_graph>
          adjacency_iterator;
  typedef void in_adjacency_iterator;
  typedef detail::csr_edge_iterator<compressed_sparse_row_graph> edge_iterator;
  typedef detail::csr_out_edge_iterator<compressed_sparse_row_graph>
          out_edge_iterator;
  typedef void in_edge_iterator;
  typedef detail::csr_thread_vertex_iterator<compressed_sparse_row_graph>
          thread_vertex_iterator;
  typedef detail::csr_thread_adjacency_iterator<compressed_sparse_row_graph>
          thread_adjacency_iterator;
  typedef void thread_in_adjacency_iterator;
  typedef detail::csr_thread_edge_iterator<compressed_sparse_row_graph>
          thread_edge_iterator;
  typedef detail::csr_thread_out_edge_iterator<compressed_sparse_row_graph>
          thread_out_edge_iterator;
  typedef void thread_in_edge_iterator;
  typedef DIRECTION directed_category;
  typedef vector_thread_iterators iterator_category;

  compressed_sparse_row_graph() : m(0), n(0), index(0),
                                  src_points(0), end_points(0),
                                  original_ids(0), internal_ids(0),
                                  initialized_by_mmap(false) {}

  compressed_sparse_row_graph(const compressed_sparse_row_graph& g)
  {
    initialized_by_mmap = false;
    deep_copy(g);
  }

  ~compressed_sparse_row_graph() { clear(); }

  compressed_sparse_row_graph& operator=(const compressed_sparse_row_graph& rhs)
  {
    clear();
    deep_copy(rhs);
    return *this;
  }

  void clear()
  {
    if (initialized_by_mmap)
    {
      initialized_by_mmap = false;
    }
    else
    {
      if (index) free(index);
      if (src_points) free(src_points);
      if (end_points) free(end_points);
      if (original_ids) free(original_ids);
      if (internal_ids) free(internal_ids);
    }

    m = 0;
    n = 0;
    index = 0;
    src_points = 0;
    end_points = 0;
    original_ids = 0;
    internal_ids = 0;
  }

  void init(size_type order, size_type size, size_type* srcs, size_type* dests)
  {
    n = order;
    m = size;
    size_type a_size = DIRECTION::is_directed() ? m : 2 * m;

    size_type* degree = (size_type*) calloc(order, sizeof(size_type));
    size_type* numEdges = (size_type*) malloc((order + 1) * sizeof(size_type));

    // Count the out degree of each vertex.
#ifdef USING_QTHREADS
    detail::increment_degrees<size_type> ids(degree, srcs);
    qt_loop_balance(0, size, ids);
#else
    #pragma mta assert parallel
    #pragma mta block schedule
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_type i = 0; i < size; ++i) mt_incr(degree[srcs[i]], 1);
#endif

    if (!DIRECTION::is_directed())
    {
#ifdef USING_QTHREADS
      detail::increment_degrees<size_type> idd(degree, dests);
      qt_loop_balance(0, size, idd);
#else
      #pragma mta assert parallel
      #pragma mta block schedule
      #ifdef _OPENMP
      #pragma omp parallel for
      #endif
      for (size_type i = 0; i < size; ++i) mt_incr(degree[dests[i]], 1);
#endif
    }

    // Find the starting index in the endpoints array for each vertex.
    numEdges[0] = 0;

    #pragma mta block schedule
    for (size_type i = 1; i <= order; ++i)
    {
      numEdges[i] = numEdges[i - 1] + degree[i - 1];
    }

    // Reset degree to be all 0's.
#ifdef USING_QTHREADS
    detail::reset_degrees<size_type> rd(degree);
    qt_loop_balance(0, order, rd);
#else
    #pragma mta block schedule
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_type i = 0; i < order; ++i) degree[i] = 0;
#endif

    size_type* srcV = (size_type*) malloc(a_size * sizeof(size_type));
    size_type* endV = (size_type*) malloc(a_size * sizeof(size_type));
    size_type* loc_original_ids =
      (size_type*) malloc(a_size * sizeof(size_type));
    size_type* loc_internal_ids = (size_type*) malloc(m * sizeof(size_type));

#ifdef USING_QTHREADS
    detail::add_csr_edges<size_type>
      ace(degree, srcs, dests, numEdges, srcV, endV,
          loc_original_ids, loc_internal_ids);
    qt_loop_balance(0, size, ace);
#else
    #pragma mta assert parallel
    #pragma mta block schedule
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_type i = 0; i < size; ++i)
    {
      size_type u = srcs[i];
      size_type vpos = numEdges[u] + mt_incr(degree[u], 1);
      srcV[vpos] = srcs[i];
      endV[vpos] = dests[i];
      loc_original_ids[vpos] = i;
      loc_internal_ids[i] = vpos;
    }
#endif

    if (!DIRECTION::is_directed())
    {
#ifdef USING_QTHREADS
      detail::add_csr_rev_edges<size_type>
        acre(degree, srcs, dests, numEdges, srcV, endV, loc_original_ids);
      qt_loop_balance(0, size, acre);
#else
      #pragma mta assert parallel
      #pragma mta block schedule
      #ifdef _OPENMP
      #pragma omp parallel for
      #endif
      for (size_type i = 0; i < size; ++i)
      {
        size_type u = dests[i];
        size_type vpos = numEdges[u] + mt_incr(degree[u], 1);
        srcV[vpos] = dests[i];
        endV[vpos] = srcs[i];
        loc_original_ids[vpos] = i;
      }
#endif
    }

    index = numEdges;
    src_points = srcV;
    end_points = endV;
    original_ids = loc_original_ids;
    internal_ids = loc_internal_ids;

    free(degree);
  }

  size_type get_order() const { return n; }
  size_type get_size() const { return m; }

  size_type get_degree(size_type i) const { return get_out_degree(i); }

  size_type get_out_degree(size_type i) const
  {
    return index[i + 1] - index[i];
  }

  vertex_iterator vertices() const { return vertex_iterator(); }

  adjacency_iterator adjacent_vertices(const vertex_descriptor& v) const
  {
    return adjacency_iterator(index, end_points, v);
  }

  edge_iterator edges() const
  {
    return edge_iterator(src_points, end_points, internal_ids);
  }

  out_edge_iterator out_edges(const vertex_descriptor& v) const
  {
    return out_edge_iterator(index, end_points, original_ids, v);
  }

  thread_vertex_iterator thread_vertices(size_type pos) const
  {
    return thread_vertex_iterator(pos);
  }

  thread_adjacency_iterator
  thread_adjacent_vertices(const vertex_descriptor& v, size_type pos) const
  {
    return thread_adjacency_iterator(index, end_points, v, pos);
  }

  thread_edge_iterator thread_edges(size_type pos) const
  {
    return thread_edge_iterator(src_points, end_points, internal_ids, pos);
  }

  thread_out_edge_iterator
  thread_out_edges(const vertex_descriptor& v, size_type pos) const
  {
    return thread_out_edge_iterator(index, end_points, original_ids, v, pos);
  }

  void print() const
  {
    printf("num Verts = %lu\n", n);
    printf("num Edges = %lu\n", m);

    for (size_type i = 0; i < n; ++i)
    {
      size_type deg = index[i + 1] - index[i];

      printf("%lu\t: [%4lu] {", i, deg);

      for (size_type j = 0; j < deg; ++j)
      {
        size_type idx = index[i];
        printf(" %lu", end_points[idx + j]);
      }

      printf(" }\n");
    }
  }

  unsigned long get_mmap_size()
  {
    // Total size:
    //
    // HEADER            TYPE              NUM VARS
    //   mmap type       unsigned long         1
    //   mmap size       unsigned long         1
    //   order           unsigned long         1
    //   size            unsigned long         1
    //
    // BODY              TYPE              DIRECTED        UNDIRECTED
    //   index           size_type           n + 1            n + 1
    //   src_points      size_type             m               2m
    //   end_points      size_type             m               2m
    //   original_ids    size_type             m               2m
    //   internal_ids    size_type             m                m

    size_type a_size = DIRECTION::is_directed() ? m : 2 * m;

    return 4 * sizeof(unsigned long) +
           ((n + 1) + 3 * a_size + m) * sizeof(size_type);
  }

  void write_mmap(void* mapped_mem)
  {
    unsigned long* ul_mapped_mem = reinterpret_cast<unsigned long*>(mapped_mem);

    size_type a_size = DIRECTION::is_directed() ? m : 2 * m;

    // Write the mmap type, mmap size, graph size, and graph order to the
    // mapped memory.
    ul_mapped_mem[0] = mmap_traits<compressed_sparse_row_graph>::type;
    ul_mapped_mem[1] = get_mmap_size();
    ul_mapped_mem[2] = n;
    ul_mapped_mem[3] = m;

    // Get pointers to the locations in the mapped memory for the arrays.
    size_type* index_ptr = reinterpret_cast<size_type*>(ul_mapped_mem + 4);
    size_type* src_points_ptr = index_ptr + (n + 1);
    size_type* end_points_ptr = src_points_ptr + a_size;
    size_type* original_ids_ptr = end_points_ptr + a_size;
    size_type* internal_ids_ptr = original_ids_ptr + a_size;

    // Write the values for the arrays to the mapped memory.
    memcpy(index_ptr, index, (n + 1) * sizeof(size_type));
    memcpy(src_points_ptr, src_points, a_size * sizeof(size_type));
    memcpy(end_points_ptr, end_points, a_size * sizeof(size_type));
    memcpy(original_ids_ptr, original_ids, a_size * sizeof(size_type));
    memcpy(internal_ids_ptr, internal_ids, m * sizeof(size_type));
  }

  void read_mmap(void* mapped_mem)
  {
    clear();

    unsigned long* ul_mapped_mem = reinterpret_cast<unsigned long*>(mapped_mem);

    // Set the size and order of the graph.
    n = ul_mapped_mem[2];
    m = ul_mapped_mem[3];

    size_type a_size = DIRECTION::is_directed() ? m : 2 * m;

    // Get pointers to the locations in the mapped memory for the arrays.
    size_type* index_ptr = reinterpret_cast<size_type*>(ul_mapped_mem + 4);
    size_type* src_points_ptr = index_ptr + (n + 1);
    size_type* end_points_ptr = src_points_ptr + a_size;
    size_type* original_ids_ptr = end_points_ptr + a_size;
    size_type* internal_ids_ptr = original_ids_ptr + a_size;

    // Set the pointers to the arrays for the graph to point to the mapped
    // memory.
    index = index_ptr;
    src_points = src_points_ptr;
    end_points = end_points_ptr;
    original_ids = original_ids_ptr;
    internal_ids = internal_ids_ptr;

    initialized_by_mmap = true;
  }

  // These accessor methods are not for normal users as they expose the
  // internals of the graph object.  They are here only because MEGRAPHS
  // needs them.
  size_type* get_index() const { return index; }
  size_type* get_src_points() const { return src_points; }
  size_type* get_end_points() const { return end_points; }
  size_type* get_original_ids() const { return original_ids; }
  size_type* get_internal_ids() const { return internal_ids; }

private:
  void deep_copy(const compressed_sparse_row_graph& rhs)
  {
    m = rhs.m;
    n = rhs.n;

    size_type a_size = DIRECTION::is_directed() ? m : 2 * m;

    index = (size_type*) malloc((n + 1) * sizeof(size_type));
    src_points = (size_type*) malloc(a_size * sizeof(size_type));
    end_points = (size_type*) malloc(a_size * sizeof(size_type));
    original_ids = (size_type*) malloc(a_size * sizeof(size_type));
    internal_ids = (size_type*) malloc(m * sizeof(size_type));

#ifdef USING_QTHREADS
    detail::arr_copy_loop<size_type> index_cpy(index, rhs.index);
    qt_loop_balance(0, n + 1, index_cpy);
    detail::arr_copy_loop<size_type> spnts_cpy(src_points, rhs.src_points);
    qt_loop_balance(0, a_size, spnts_cpy);
    detail::arr_copy_loop<size_type> epnts_cpy(end_points, rhs.end_points);
    qt_loop_balance(0, a_size, epnts_cpy);
    detail::arr_copy_loop<size_type> oids_cpy(original_ids, rhs.original_ids);
    qt_loop_balance(0, a_size, oids_cpy);
    detail::arr_copy_loop<size_type> iids_cpy(internal_ids, rhs.internal_ids);
    qt_loop_balance(0, m, iids_cpy);
#else
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_type i = 0; i < n + 1; ++i) index[i] = rhs.index[i];

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_type i = 0; i < a_size; ++i) src_points[i] = rhs.src_points[i];

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_type i = 0; i < a_size; ++i) end_points[i] = rhs.end_points[i];

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_type i = 0; i < a_size; ++i)
    {
      original_ids[i] = rhs.original_ids[i];
    }

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_type i = 0; i < m; ++i) internal_ids[i] = rhs.internal_ids[i];
#endif
  }

private:
  size_type m;      // # of edges
  size_type n;      // # of vertices

  // Indexes into the edge-arrays for the adjacencies.  Its size is |V|+1.  A
  // vertex, v, has neighbors in {src,end}_points from {src,end}_points[v] to
  // {src,end}_points[v+1].
  size_type* index;

  // src_points and end_points store the source and target vertices for each
  // edge, respectively.  These arrays are of size |E| if directed and 2*|E|
  // if undirected.
  size_type* src_points;
  size_type* end_points;

  // The order of the edges is rearranged from the passed in order during the
  // creation of the graph.  However, the ids of the edges should be assigned
  // based on the original order not the new order.  We need a way to map the
  // original order of the edges to the internal order and vice versa in order
  // to map the correct ids to the correct edges.
  size_type* original_ids;
  size_type* internal_ids;

  bool initialized_by_mmap;
};

/***/

template <typename DIRECTION>
inline
typename compressed_sparse_row_graph<DIRECTION>::size_type
num_vertices(const compressed_sparse_row_graph<DIRECTION>& g)
{
  return g.get_order();
}

/***/

template <typename DIRECTION>
inline
typename compressed_sparse_row_graph<DIRECTION>::size_type
num_edges(const compressed_sparse_row_graph<DIRECTION>& g)
{
  return g.get_size();
}

/***/

template <typename DIRECTION>
inline
typename compressed_sparse_row_graph<DIRECTION>::vertex_descriptor
source(
    const typename compressed_sparse_row_graph<DIRECTION>::edge_descriptor& e,
    const compressed_sparse_row_graph<DIRECTION>& g)
{
  return e.first;
}

/***/

template <typename DIRECTION>
inline
typename compressed_sparse_row_graph<DIRECTION>::vertex_descriptor
target(
    const typename compressed_sparse_row_graph<DIRECTION>::edge_descriptor& e,
    const compressed_sparse_row_graph<DIRECTION>& g)
{
  return e.second;
}

/***/

template <typename DIRECTION>
inline
typename compressed_sparse_row_graph<DIRECTION>::size_type
degree(
    const typename compressed_sparse_row_graph<DIRECTION>::vertex_descriptor& v,
    const compressed_sparse_row_graph<DIRECTION>& g)
{
  return g.get_degree(v);
}

/***/

template <typename DIRECTION>
inline
typename compressed_sparse_row_graph<DIRECTION>::size_type
out_degree(
    const typename compressed_sparse_row_graph<DIRECTION>::vertex_descriptor& v,
    const compressed_sparse_row_graph<DIRECTION>& g)
{
  return g.get_out_degree(v);
}

/***/

template <typename DIRECTION>
inline
typename compressed_sparse_row_graph<DIRECTION>::vertex_iterator
vertices(const compressed_sparse_row_graph<DIRECTION>& g)
{
  return g.vertices();
}

/***/

template <typename DIRECTION>
inline
typename compressed_sparse_row_graph<DIRECTION>::adjacency_iterator
adjacent_vertices(
    const typename compressed_sparse_row_graph<DIRECTION>::vertex_descriptor& v,
    const compressed_sparse_row_graph<DIRECTION>& g)
{
  return g.adjacent_vertices(v);
}

/***/

template <typename DIRECTION>
inline
typename compressed_sparse_row_graph<DIRECTION>::edge_iterator
edges(const compressed_sparse_row_graph<DIRECTION>& g)
{
  return g.edges();
}

/***/

template <typename DIRECTION>
inline
typename compressed_sparse_row_graph<DIRECTION>::out_edge_iterator
out_edges(
    const typename compressed_sparse_row_graph<DIRECTION>::vertex_descriptor& v,
    const compressed_sparse_row_graph<DIRECTION>& g)
{
  return g.out_edges(v);
}

/***/

template <typename DIRECTION>
inline
typename compressed_sparse_row_graph<DIRECTION>::thread_vertex_iterator
thread_vertices(typename compressed_sparse_row_graph<DIRECTION>::size_type pos,
                const compressed_sparse_row_graph<DIRECTION>& g)
{
  return g.thread_vertices(pos);
}

/***/

template <typename DIRECTION>
inline
typename compressed_sparse_row_graph<DIRECTION>::thread_adjacency_iterator
thread_adjacent_vertices(
  const typename compressed_sparse_row_graph<DIRECTION>::vertex_descriptor& v,
  typename compressed_sparse_row_graph<DIRECTION>::size_type pos,
  const compressed_sparse_row_graph<DIRECTION>& g)
{
  return g.thread_adjacent_vertices(v, pos);
}

/***/

template <typename DIRECTION>
inline
typename compressed_sparse_row_graph<DIRECTION>::thread_edge_iterator
thread_edges(typename compressed_sparse_row_graph<DIRECTION>::size_type pos,
             const compressed_sparse_row_graph<DIRECTION>& g)
{
  return g.thread_edges(pos);
}

/***/

template <typename DIRECTION>
inline
typename compressed_sparse_row_graph<DIRECTION>::thread_out_edge_iterator
thread_out_edges(
  const typename compressed_sparse_row_graph<DIRECTION>::vertex_descriptor& v,
  typename compressed_sparse_row_graph<DIRECTION>::size_type pos,
  const compressed_sparse_row_graph<DIRECTION>& g)
{
  return g.thread_out_edges(v, pos);
}

/***/

template <typename DIRECTION>
inline
typename compressed_sparse_row_graph<DIRECTION>::vertex_descriptor
null_vertex(const compressed_sparse_row_graph<DIRECTION>& g)
{
  return (std::numeric_limits<typename compressed_sparse_row_graph<DIRECTION>::
                                       vertex_descriptor>::max)();
}

/***/

template <typename DIRECTION>
inline
typename compressed_sparse_row_graph<DIRECTION>::edge_descriptor
null_edge(const compressed_sparse_row_graph<DIRECTION>& g)
{
  return typename compressed_sparse_row_graph<DIRECTION>::edge_descriptor();
}

/***/

template <typename ITERATOR, typename DIRECTION>
inline
bool
is_valid(ITERATOR& iter,
         typename compressed_sparse_row_graph<DIRECTION>::size_type p,
         const compressed_sparse_row_graph<DIRECTION>& tg)
{
  return true;
}

/***/

template <typename DIRECTION>
inline
bool
is_directed(const compressed_sparse_row_graph<DIRECTION>& g)
{
  return DIRECTION::is_directed();
}

/***/

template <typename DIRECTION>
inline
bool
is_undirected(const compressed_sparse_row_graph<DIRECTION>& g)
{
  return !is_directed(g);
}

/***/

template <typename DIRECTION>
inline
bool
is_bidirectional(const compressed_sparse_row_graph<DIRECTION>& g)
{
  return DIRECTION::is_bidirectional();
}

/***/

template <typename DIRECTION>
class vertex_id_map<compressed_sparse_row_graph<DIRECTION> > :
  public put_get_helper<
             typename compressed_sparse_row_graph<DIRECTION>::size_type,
             vertex_id_map<compressed_sparse_row_graph<DIRECTION> > > {
public:
  typedef typename compressed_sparse_row_graph<DIRECTION>::vertex_descriptor
          key_type;
  typedef typename compressed_sparse_row_graph<DIRECTION>::size_type value_type;

  vertex_id_map() {}
  value_type operator[] (const key_type& k) const { return k; }
};


template <typename DIRECTION>
class vertex_id_map<const compressed_sparse_row_graph<DIRECTION> > :
  public put_get_helper<
             typename compressed_sparse_row_graph<DIRECTION>::size_type,
             vertex_id_map<const compressed_sparse_row_graph<DIRECTION> > > {
public:
  typedef typename compressed_sparse_row_graph<DIRECTION>::vertex_descriptor
          key_type;
  typedef typename compressed_sparse_row_graph<DIRECTION>::size_type value_type;

  vertex_id_map() {}
  value_type operator[] (const key_type& k) const { return k; }
};

/***/

template <typename DIRECTION>
class edge_id_map<compressed_sparse_row_graph<DIRECTION> > :
  public put_get_helper<
             typename compressed_sparse_row_graph<DIRECTION>::size_type,
             edge_id_map<compressed_sparse_row_graph<DIRECTION> > > {
public:
  typedef typename compressed_sparse_row_graph<DIRECTION>::edge_descriptor
          key_type;
  typedef typename compressed_sparse_row_graph<DIRECTION>::size_type value_type;

  edge_id_map() {}
  value_type operator[] (const key_type& k) const { return k.id; }
};

template <typename DIRECTION>
class edge_id_map<const compressed_sparse_row_graph<DIRECTION> > :
  public put_get_helper<
             typename compressed_sparse_row_graph<DIRECTION>::size_type,
             edge_id_map<const compressed_sparse_row_graph<DIRECTION> > > {
public:
  typedef typename compressed_sparse_row_graph<DIRECTION>::edge_descriptor
          key_type;
  typedef typename compressed_sparse_row_graph<DIRECTION>::size_type value_type;

  edge_id_map() {}
  value_type operator[] (const key_type& k) const { return k.id; }
};

/***/

template <typename Graph>
struct default_hash_func<detail::csr_edge_adapter<Graph> > {
  hash_size_type operator()(const detail::csr_edge_adapter<Graph>& key) const
  { return integer_hash_func(key.id); }
};

/***/

template <typename DIRECTION>
class mmap_traits<compressed_sparse_row_graph<DIRECTION> > {
public:
  static const unsigned long type = 100 + DIRECTION::direction_type;
};

/***/

template <typename DIRECTION>
inline
void init(typename compressed_sparse_row_graph<DIRECTION>::size_type n,
          typename compressed_sparse_row_graph<DIRECTION>::size_type m,
          typename compressed_sparse_row_graph<DIRECTION>::size_type* srcs,
          typename compressed_sparse_row_graph<DIRECTION>::size_type* dests,
          compressed_sparse_row_graph<DIRECTION>& g)
{
  return g.init(n, m, srcs, dests);
}

/***/

template <typename DIRECTION>
inline
void clear(compressed_sparse_row_graph<DIRECTION>& g)
{
  return g.clear();
}

/***************************************/

namespace detail {

template <typename Graph>
class csr_in_edge_iterator {
private:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

public:
  csr_in_edge_iterator() : adj(0), original_ids(0), vid(0) {}

  csr_in_edge_iterator(size_type* ind, size_type* adj_, size_type* oi,
                       size_type i) :
    adj(adj_ + ind[i]), original_ids(oi + ind[i]), vid(i) {}

  edge_descriptor operator[](size_type p) const
  {
    return edge_descriptor(adj[p], vid, original_ids[p]);
  }

private:
  size_type* adj;
  size_type* original_ids;
  size_type vid;
};

/***/

template <typename Graph>
class csr_thread_in_edge_iterator {
private:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

public:
  csr_thread_in_edge_iterator() : adj(0), original_ids(0), vid(0) {}

  csr_thread_in_edge_iterator(size_type* ind, vertex_descriptor* adj_,
                               size_type* oi, size_type i, size_type adjp) :
    adj(adj_ + ind[i] + adjp), original_ids(oi + ind[i] + adjp), vid(i) {}

  csr_thread_in_edge_iterator(const csr_thread_in_edge_iterator& rhs) :
    adj(rhs.adj), original_ids(rhs.original_ids), vid(rhs.vid) {}

  csr_thread_in_edge_iterator&
  operator=(const csr_thread_in_edge_iterator& rhs)
  {
    adj = rhs.adj;
    original_ids = rhs.original_ids;
    vid = rhs.vid;
    return *this;
  }

  csr_thread_in_edge_iterator& operator++()
  {
    ++adj;
    ++original_ids;
    return *this;
  }

  csr_thread_in_edge_iterator& operator++(int)
  {
    csr_thread_in_edge_iterator temp(*this);

    ++adj;
    ++original_ids;
    return temp;
  }

  edge_descriptor operator*() const
  {
    return edge_descriptor(*adj, vid, *original_ids);
  }

  bool operator==(const csr_thread_in_edge_iterator& rhs)
  { return adj == rhs.adj; }
  bool operator!=(const csr_thread_in_edge_iterator& rhs)
  { return adj != rhs.adj; }

private:
  vertex_descriptor* adj;
  size_type* original_ids;
  size_type vid;
};

}

/***/

template <>
class compressed_sparse_row_graph<bidirectionalS> {
public:
  typedef unsigned long size_type;
  typedef unsigned long vertex_descriptor;
  typedef detail::csr_edge_adapter<compressed_sparse_row_graph> edge_descriptor;
  typedef detail::csr_vertex_iterator<compressed_sparse_row_graph>
          vertex_iterator;
  typedef detail::csr_adjacency_iterator<compressed_sparse_row_graph>
          adjacency_iterator;
  typedef detail::csr_adjacency_iterator<compressed_sparse_row_graph>
          in_adjacency_iterator;
  typedef detail::csr_edge_iterator<compressed_sparse_row_graph> edge_iterator;
  typedef detail::csr_out_edge_iterator<compressed_sparse_row_graph>
          out_edge_iterator;
  typedef detail::csr_in_edge_iterator<compressed_sparse_row_graph>
          in_edge_iterator;
  typedef detail::csr_thread_vertex_iterator<compressed_sparse_row_graph>
          thread_vertex_iterator;
  typedef detail::csr_thread_adjacency_iterator<compressed_sparse_row_graph>
          thread_adjacency_iterator;
  typedef detail::csr_thread_adjacency_iterator<compressed_sparse_row_graph>
          thread_in_adjacency_iterator;
  typedef detail::csr_thread_edge_iterator<compressed_sparse_row_graph>
          thread_edge_iterator;
  typedef detail::csr_thread_out_edge_iterator<compressed_sparse_row_graph>
          thread_out_edge_iterator;
  typedef detail::csr_thread_in_edge_iterator<compressed_sparse_row_graph>
          thread_in_edge_iterator;
  typedef bidirectionalS directed_category;
  typedef vector_thread_iterators iterator_category;

  compressed_sparse_row_graph() : m(0), n(0), index(0), rev_index(0),
                                  src_points(0), rev_src_points(0),
                                  end_points(0), rev_end_points(0),
                                  original_ids(0), rev_original_ids(0),
                                  internal_ids(0), initialized_by_mmap(false) {}

  compressed_sparse_row_graph(const compressed_sparse_row_graph& g)
  {
    initialized_by_mmap = false;
    deep_copy(g);
  }

  ~compressed_sparse_row_graph() { clear(); }

  compressed_sparse_row_graph& operator=(const compressed_sparse_row_graph& rhs)
  {
    clear();
    deep_copy(rhs);
    return *this;
  }

  void clear()
  {
    if (initialized_by_mmap)
    {
      initialized_by_mmap = false;
    }
    else
    {
      if (index) free(index);
      if (rev_index) free(rev_index);
      if (src_points) free(src_points);
      if (rev_src_points) free(rev_src_points);
      if (end_points) free(end_points);
      if (rev_end_points) free(rev_end_points);
      if (original_ids) free(original_ids);
      if (rev_original_ids) free(rev_original_ids);
      if (internal_ids) free(internal_ids);
    }

    m = 0;
    n = 0;
    index = 0;
    rev_index = 0;
    src_points = 0;
    rev_src_points = 0;
    end_points = 0;
    rev_end_points = 0;
    original_ids = 0;
    rev_original_ids = 0;
    internal_ids = 0;
  }

  void init(size_type order, size_type size, size_type* srcs, size_type* dests)
  {
    clear();

    n = order;
    m = size;
    size_type a_size = m;

    size_type* out_degree = (size_type*) calloc(order, sizeof(size_type));
    size_type* in_degree = (size_type*) calloc(order, sizeof(size_type));
    size_type* out_numEdges =
      (size_type*) malloc((order + 1) * sizeof(size_type));
    size_type* in_numEdges =
      (size_type*) malloc((order + 1) * sizeof(size_type));

#ifdef USING_QTHREADS
    detail::increment_degrees<size_type> ids(out_degree, srcs);
    qt_loop_balance(0, size, ids);

    detail::increment_degrees<size_type> idd(in_degree, dests);
    qt_loop_balance(0, size, idd);
#else
    #pragma mta assert parallel
    #pragma mta block schedule
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_type i = 0; i < size; ++i)
    {
      mt_incr(out_degree[srcs[i]], 1);
      mt_incr(in_degree[dests[i]], 1);
    }
#endif

    in_numEdges[0] = 0;
    out_numEdges[0] = 0;

    #pragma mta block schedule
    for (size_type i = 1; i <= order; ++i)
    {
      out_numEdges[i] = out_numEdges[i - 1] + out_degree[i - 1];
    }

    #pragma mta block schedule
    for (size_type i = 1; i <= order; ++i)
    {
      in_numEdges[i] = in_numEdges[i - 1] + in_degree[i - 1];
    }

    // Reset out_degree and in_degree to be all 0's.
#ifdef USING_QTHREADS
    detail::reset_degrees<size_type> rdo(out_degree);
    qt_loop_balance(0, order, rdo);

    detail::reset_degrees<size_type> rdi(in_degree);
    qt_loop_balance(0, order, rdi);
#else
    #pragma mta block schedule
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_type i = 0; i < order; ++i)
    {
      out_degree[i] = 0;
      in_degree[i] = 0;
    }
#endif

    size_type* out_endV = (size_type*) malloc(a_size * sizeof(size_type));
    size_type* out_srcV = (size_type*) malloc(a_size * sizeof(size_type));
    size_type* in_endV = (size_type*) malloc(a_size * sizeof(size_type));
    size_type* in_srcV = (size_type*) malloc(a_size * sizeof(size_type));
    size_type* cross_cross_index =
      (size_type*) malloc(a_size * sizeof(size_type));
    size_type* loc_original_ids =
      (size_type*) malloc(a_size * sizeof(size_type));
    size_type* loc_internal_ids = (size_type*) malloc(m * sizeof(size_type));

#ifdef USING_QTHREADS
    detail::add_csr_bid_edges<size_type>
      acbe(out_degree, in_degree, srcs, dests, out_numEdges, in_numEdges,
           out_srcV, out_endV, in_srcV, in_endV,
           loc_original_ids, loc_internal_ids, cross_cross_index);
    qt_loop_balance(0, size, acbe);
#else
    #pragma mta assert parallel
    #pragma mta block schedule
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_type i = 0; i < size; ++i)
    {
      size_type src = srcs[i];
      size_type dest = dests[i];

      // Add the edge to the source vertex.
      size_type out_vp = out_numEdges[src] + mt_incr(out_degree[src], 1);
      out_endV[out_vp] = dest;
      out_srcV[out_vp] = src;
      loc_original_ids[out_vp] = i;
      loc_internal_ids[i] = out_vp;

      // Add the edge to the dest vertex.
      size_type in_vp = in_numEdges[dest] + mt_incr(in_degree[dest], 1);
      in_endV[in_vp] = src;
      in_srcV[in_vp] = dest;
      cross_cross_index[in_vp] = out_vp;
    }
#endif

    // Set the mapping from the reverse edges to the original ids.
    size_type* loc_rev_original_ids =
      (size_type*) malloc(a_size * sizeof(size_type));

#ifdef USING_QTHREADS
    detail::add_rev_original_ids<size_type>
      aroi(loc_original_ids, loc_rev_original_ids, cross_cross_index);
    qt_loop_balance(0, a_size, aroi);
#else
    #pragma mta assert nodep
    #pragma mta block schedule
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_type i = 0; i < a_size; ++i)
    {
      loc_rev_original_ids[i] = loc_original_ids[cross_cross_index[i]];
    }
#endif

    index = out_numEdges;
    rev_index = in_numEdges;
    src_points = out_srcV;
    rev_src_points = in_srcV;
    end_points = out_endV;
    rev_end_points = in_endV;
    original_ids = loc_original_ids;
    rev_original_ids = loc_rev_original_ids;
    internal_ids = loc_internal_ids;

    free(out_degree);
    free(in_degree);
    free(cross_cross_index);
  }

  size_type get_order() const { return n; }
  size_type get_size() const { return m; }

  size_type get_degree(vertex_descriptor i) const
  {
    return get_in_degree(i) + get_out_degree(i);
  }

  size_type get_out_degree(vertex_descriptor i) const
  {
    return index[i + 1] - index[i];
  }

  size_type get_in_degree(vertex_descriptor i) const
  {
    return rev_index[i + 1] - rev_index[i];
  }

  vertex_iterator vertices() const { return vertex_iterator(); }

  adjacency_iterator adjacent_vertices(const vertex_descriptor& v) const
  {
    return adjacency_iterator(index, end_points, v);
  }

  in_adjacency_iterator in_adjacent_vertices(const vertex_descriptor& v) const
  {
    return in_adjacency_iterator(rev_index, rev_end_points, v);
  }

  edge_iterator edges() const
  {
    return edge_iterator(src_points, end_points, internal_ids);
  }

  out_edge_iterator out_edges(const vertex_descriptor& v) const
  {
    return out_edge_iterator(index, end_points, original_ids, v);
  }

  in_edge_iterator in_edges(const vertex_descriptor& v) const
  {
    return in_edge_iterator(rev_index, rev_end_points, rev_original_ids, v);
  }

  thread_vertex_iterator thread_vertices(size_type pos) const
  {
    return thread_vertex_iterator(pos);
  }

  thread_adjacency_iterator
  thread_adjacent_vertices(const vertex_descriptor& v, size_type pos) const
  {
    return thread_adjacency_iterator(index, end_points, v, pos);
  }

  thread_in_adjacency_iterator
  thread_in_adjacent_vertices(const vertex_descriptor& v, size_type pos) const
  {
    return thread_in_adjacency_iterator(rev_index, rev_end_points, v, pos);
  }

  thread_edge_iterator thread_edges(size_type pos) const
  {
    return thread_edge_iterator(src_points, end_points, internal_ids, pos);
  }

  thread_out_edge_iterator
  thread_out_edges(const vertex_descriptor& v, size_type pos) const
  {
    return thread_out_edge_iterator(index, end_points, original_ids, v, pos);
  }

  thread_in_edge_iterator
  thread_in_edges(const vertex_descriptor& v, size_type pos) const
  {
    return thread_in_edge_iterator(rev_index, rev_end_points, rev_original_ids,
                                   v, pos);
  }

  void print() const
  {
    printf("num Verts = %lu\n", n);
    printf("num Edges = %lu\n", m);

    for (size_type i = 0; i < n; ++i)
    {
      size_type deg = index[i + 1] - index[i];
      printf("%lu\t: [%4lu] {", i, deg);

      for (size_type j = 0; j < deg; ++j)
      {
        size_type idx = index[i];
        printf(" %lu", end_points[idx + j]);
      }

      printf(" }\n");
    }
  }

  unsigned long get_mmap_size()
  {
    // Total size:
    //
    // HEADER                TYPE               NUM VARS
    //   mmap type           unsigned long          1
    //   mmap size           unsigned long          1
    //   order               unsigned long          1
    //   size                unsigned long          1
    //
    // BODY                  TYPE             BIDIRECTIONAL
    //   index               size_type            n + 1
    //   rev_index           size_type            n + 1
    //   src_points          size_type              m
    //   rev_src_points      size_type              m
    //   end_points          size_type              m
    //   rev_end_points      size_type              m
    //   original_ids        size_type              m
    //   rev_original_ids    size_type              m
    //   internal_ids        size_type              m

    return 4 * sizeof(unsigned long) +
           (2 * (n + 1) + 7 * m) * sizeof(size_type);
  }

  void write_mmap(void* mapped_mem)
  {
    unsigned long* ul_mapped_mem = reinterpret_cast<unsigned long*>(mapped_mem);

    // Write the mmap type, mmap size, graph size, and graph order to the
    // mapped memory.
    ul_mapped_mem[0] = mmap_traits<compressed_sparse_row_graph>::type;
    ul_mapped_mem[1] = get_mmap_size();
    ul_mapped_mem[2] = n;
    ul_mapped_mem[3] = m;

    // Get pointers to the locations in the mapped memory for the arrays.
    size_type* index_ptr = reinterpret_cast<size_type*>(ul_mapped_mem + 4);
    size_type* rev_index_ptr = index_ptr + (n + 1);
    size_type* src_points_ptr = rev_index_ptr + (n + 1);
    size_type* rev_src_points_ptr = src_points_ptr + m;
    size_type* end_points_ptr = rev_src_points_ptr + m;
    size_type* rev_end_points_ptr = end_points_ptr + m;
    size_type* original_ids_ptr = rev_end_points_ptr + m;
    size_type* rev_original_ids_ptr = original_ids_ptr + m;
    size_type* internal_ids_ptr = rev_original_ids_ptr + m;

    // Write the values for the arrays to the mapped memory.
    memcpy(index_ptr, index, (n + 1) * sizeof(size_type));
    memcpy(rev_index_ptr, rev_index, (n + 1) * sizeof(size_type));
    memcpy(src_points_ptr, src_points, m * sizeof(size_type));
    memcpy(rev_src_points_ptr, rev_src_points, m * sizeof(size_type));
    memcpy(end_points_ptr, end_points, m * sizeof(size_type));
    memcpy(rev_end_points_ptr, rev_end_points, m * sizeof(size_type));
    memcpy(original_ids_ptr, original_ids, m * sizeof(size_type));
    memcpy(rev_original_ids_ptr, rev_original_ids, m * sizeof(size_type));
    memcpy(internal_ids_ptr, internal_ids, m * sizeof(size_type));
  }

  void read_mmap(void* mapped_mem)
  {
    clear();

    unsigned long* ul_mapped_mem = reinterpret_cast<unsigned long*>(mapped_mem);

    // Set the size and order of the graph.
    n = ul_mapped_mem[2];
    m = ul_mapped_mem[3];

    // Get pointers to the locations in the mapped memory for the arrays.
    size_type* index_ptr = reinterpret_cast<size_type*>(ul_mapped_mem + 4);
    size_type* rev_index_ptr = index_ptr + (n + 1);
    size_type* src_points_ptr = rev_index_ptr + (n + 1);
    size_type* rev_src_points_ptr = src_points_ptr + m;
    size_type* end_points_ptr = rev_src_points_ptr + m;
    size_type* rev_end_points_ptr = end_points_ptr + m;
    size_type* original_ids_ptr = rev_end_points_ptr + m;
    size_type* rev_original_ids_ptr = original_ids_ptr + m;
    size_type* internal_ids_ptr = rev_original_ids_ptr + m;

    // Set the pointers to the arrays for the graph to point to the mapped
    // memory.
    index = index_ptr;
    rev_index = rev_index_ptr;
    src_points = src_points_ptr;
    rev_src_points = rev_src_points_ptr;
    end_points = end_points_ptr;
    rev_end_points = rev_end_points_ptr;
    original_ids = original_ids_ptr;
    rev_original_ids = rev_original_ids_ptr;
    internal_ids = internal_ids_ptr;

    initialized_by_mmap = true;
  }

  // These accessor methods are not for normal users as they expose the
  // internals of the graph object.  They are here only because MEGRAPHS
  // needs them.
  size_type* get_index() const { return index; }
  size_type* get_rev_index() const { return rev_index; }
  size_type* get_src_points() const { return src_points; }
  size_type* get_rev_src_points() { return rev_src_points; }
  size_type* get_end_points() const { return end_points; }
  size_type* get_rev_end_points() { return rev_end_points; }
  size_type* get_original_ids() const { return original_ids; }
  size_type* get_rev_original_ids() const { return rev_original_ids; }
  size_type* get_internal_ids() const { return internal_ids; }

private:
  void deep_copy(const compressed_sparse_row_graph& rhs)
  {
    m = rhs.m;
    n = rhs.n;

    size_type a_size = m;

    index = (size_type*) malloc((n + 1) * sizeof(size_type));
    rev_index = (size_type*) malloc((n + 1) * sizeof(size_type));
    src_points = (size_type*) malloc(a_size * sizeof(size_type));
    rev_src_points = (size_type*) malloc(a_size * sizeof(size_type));
    end_points = (size_type*) malloc(a_size * sizeof(size_type));
    rev_end_points = (size_type*) malloc(a_size * sizeof(size_type));
    original_ids = (size_type*) malloc(a_size * sizeof(size_type));
    rev_original_ids = (size_type*) malloc(a_size * sizeof(size_type));
    internal_ids = (size_type*) malloc(m * sizeof(size_type));

#ifdef USING_QTHREADS
    detail::arr_copy_loop<size_type> index_cpy(index, rhs.index);
    qt_loop_balance(0, n + 1, index_cpy);
    detail::arr_copy_loop<size_type> rindex_cpy(rev_index, rhs.rev_index);
    qt_loop_balance(0, n + 1, rindex_cpy);
    detail::arr_copy_loop<size_type> spnts_cpy(src_points, rhs.src_points);
    qt_loop_balance(0, a_size, spnts_cpy);
    detail::arr_copy_loop<size_type> rspnts_cpy(rev_src_points,
                                                rhs.rev_src_points);
    qt_loop_balance(0, a_size, rspnts_cpy);
    detail::arr_copy_loop<size_type> epnts_cpy(end_points, rhs.end_points);
    qt_loop_balance(0, a_size, epnts_cpy);
    detail::arr_copy_loop<size_type> repnts_cpy(rev_end_points,
                                                rhs.rev_end_points);
    qt_loop_balance(0, a_size, repnts_cpy);
    detail::arr_copy_loop<size_type> oids_cpy(original_ids, rhs.original_ids);
    qt_loop_balance(0, a_size, oids_cpy);
    detail::arr_copy_loop<size_type> roids_cpy(rev_original_ids,
                                               rhs.rev_original_ids);
    qt_loop_balance(0, a_size, roids_cpy);
    detail::arr_copy_loop<size_type> iids_cpy(internal_ids, rhs.internal_ids);
    qt_loop_balance(0, a_size, iids_cpy);
#else
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_type i = 0; i < n + 1; ++i) index[i] = rhs.index[i];

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_type i = 0; i < n + 1; ++i) rev_index[i] = rhs.rev_index[i];

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_type i = 0; i < a_size; ++i) src_points[i] = rhs.src_points[i];

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_type i = 0; i < a_size; ++i)
    {
      rev_src_points[i] = rhs.rev_src_points[i];
    }

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_type i = 0; i < a_size; ++i) end_points[i] = rhs.end_points[i];

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_type i = 0; i < a_size; ++i)
    {
      rev_end_points[i] = rhs.rev_end_points[i];
    }

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_type i = 0; i < a_size; ++i)
    {
      original_ids[i] = rhs.original_ids[i];
    }

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_type i = 0; i < a_size; ++i)
    {
      rev_original_ids[i] = rhs.rev_original_ids[i];
    }

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_type i = 0; i < a_size; ++i)
    {
      internal_ids[i] = rhs.internal_ids[i];
    }
#endif
  }

private:
  size_type m;  // # of edges
  size_type n;  // # of vertices

  // Indexes into the edge-arrays for the adjacencies.  Its size is |V|+1.  A
  // vertex, v, has neighbors in {src,end}_points from {src,end}_points[v] to
  // {src,end}_points[v+1].  The same is true for rev_index with respect to
  // rev_{src,end}_points.
  size_type* index;
  size_type* rev_index;

  // src_points and end_points store the source and target vertices for each
  // edge, respectively.  rev_src_points and rev_end_points store the source
  // and target vertices for each inedge, respectively.  These arrays are of
  // size |E|.
  size_type* src_points;
  size_type* rev_src_points;
  size_type* end_points;
  size_type* rev_end_points;

  // The order of the edges is rearranged from the passed in order during the
  // creation of the graph.  However, the ids of the edges should be assigned
  // based on the original order not the new order.  We need a way to map the
  // original order of the edges to the internal order (internal_ids) and
  // vice versa (original_ids) in order to map the correct ids to the
  // correct edges.  We also need a mapping from the internal order of the
  // reverse edges to the original order (rev_original_ids).
  size_type* original_ids;
  size_type* rev_original_ids;
  size_type* internal_ids;

  bool initialized_by_mmap;
};

/***/

inline
compressed_sparse_row_graph<bidirectionalS>::size_type
in_degree(
    const compressed_sparse_row_graph<bidirectionalS>::vertex_descriptor& v,
    const compressed_sparse_row_graph<bidirectionalS>& g)
{
  return g.get_in_degree(v);
}

/***/

inline
compressed_sparse_row_graph<bidirectionalS>::in_adjacency_iterator
in_adjacent_vertices(
    const compressed_sparse_row_graph<bidirectionalS>::vertex_descriptor& v,
    const compressed_sparse_row_graph<bidirectionalS>& g)
{
  return g.in_adjacent_vertices(v);
}

/***/

inline
compressed_sparse_row_graph<bidirectionalS>::in_edge_iterator
in_edges(
    const compressed_sparse_row_graph<bidirectionalS>::vertex_descriptor& v,
    const compressed_sparse_row_graph<bidirectionalS>& g)
{
  return g.in_edges(v);
}

/***/

inline
compressed_sparse_row_graph<bidirectionalS>::thread_in_adjacency_iterator
thread_in_adjacent_vertices(
    const compressed_sparse_row_graph<bidirectionalS>::vertex_descriptor& v,
    compressed_sparse_row_graph<bidirectionalS>::size_type pos,
    const compressed_sparse_row_graph<bidirectionalS>& g)
{
  return g.thread_in_adjacent_vertices(v, pos);
}

/***/

inline
compressed_sparse_row_graph<bidirectionalS>::thread_in_edge_iterator
thread_in_edges(
    const compressed_sparse_row_graph<bidirectionalS>::vertex_descriptor& v,
    compressed_sparse_row_graph<bidirectionalS>::size_type pos,
    const compressed_sparse_row_graph<bidirectionalS>& g)
{
  return g.thread_in_edges(v, pos);
}

/***/

// Needed to implement to go over both the in and out going edges if we want
// UNDIRECTED.
template <typename visitor>
void
visit_edges(compressed_sparse_row_graph<bidirectionalS>& g,
            compressed_sparse_row_graph<bidirectionalS>::vertex_descriptor v,
            visitor& vis, int par_cutoff, bool use_future, int directed)
{
  typedef compressed_sparse_row_graph<bidirectionalS>::size_type size_type;
  typedef compressed_sparse_row_graph<bidirectionalS>::edge_descriptor
          edge_descriptor;
  typedef compressed_sparse_row_graph<bidirectionalS>::out_edge_iterator
          out_edge_iterator;
  typedef compressed_sparse_row_graph<bidirectionalS>::in_edge_iterator
          in_edge_iterator;

  size_type deg;

  if (directed != REVERSED)
  {
    out_edge_iterator eit = out_edges(v, g);
    deg = out_degree(v, g);

    if (deg >= par_cutoff)
    {
      if (use_future)
      {
        #pragma mta assert parallel
        #pragma mta loop future
        for (size_type i = 0; i < deg; ++i)
        {
          edge_descriptor e = eit[i];
          vis(e);
        }
      }
      else
      {
        #pragma mta assert parallel
        for (size_type i = 0; i < deg; ++i)
        {
          edge_descriptor e = eit[i];
          vis(e);
        }
      }
    }
    else
    {
      for (size_type i = 0; i < deg; ++i)
      {
        edge_descriptor e = eit[i];
        vis(e);
      }
    }
  }

  if (directed != DIRECTED)
  {
    in_edge_iterator eit = in_edges(v, g);
    deg = in_degree(v, g);

    if (deg >= par_cutoff)
    {
      if (use_future)
      {
        #pragma mta assert parallel
        #pragma mta loop future
        for (size_type i = 0; i < deg; ++i)
        {
          edge_descriptor e = eit[i];
          vis.visit_in(e);
        }
      }
      else
      {
        #pragma mta assert parallel
        for (size_type i = 0; i < deg; ++i)
        {
          edge_descriptor e = eit[i];
          vis.visit_in(e);
        }
      }
    }
    else
    {
      for (size_type i = 0; i < deg; ++i)
      {
        edge_descriptor e = eit[i];
        vis.visit_in(e);
      }
    }
  }
}

}

#endif
