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
/*! \file neighborhoods.hpp

    \author Jon Berry (jberry@sandia.gov)

    \date 9/2008
*/
/****************************************************************************/

#ifndef MTGL_NEIGHBORHOODS_HPP
#define MTGL_NEIGHBORHOODS_HPP

#include <cstdio>
#include <climits>
#include <cstdlib>
#include <list>

#ifndef _WIN32
#include <sys/time.h>
#endif

#include <mtgl/mtgl_adapter.hpp>
#include <mtgl/util.hpp>
#include <mtgl/triangles.hpp>
#include <mtgl/ufl.h>
#include <mtgl/metrics.hpp>
#include <mtgl/xmt_hash_table.hpp>
#include <mtgl/visit_adj.hpp>

#ifdef __MTA__
#include <sys/mta_task.h>
#include <machine/runtime.h>
#endif

#define THRESH 400        // LFR settings for "Tolerating..." paper
#define FE_FACTOR 35

//#define THRESH 50
//#define FE_FACTOR 6

inline double getWtime()
{
  // Returns a double value indicating seconds elapsed since some arbitrary
  // point in the past.
#ifndef _WIN32
  timeval theTimeVal;
  gettimeofday(&theTimeVal, 0);

  return static_cast<double>(theTimeVal.tv_sec) +
         (static_cast<double>(theTimeVal.tv_usec) / 1000000);
#else
  return(0.0);
#endif
}

namespace mtgl {

template <typename T>
inline static std::list<std::list<T> > choose(T* a, int sz, int k)
{
  std::list<std::list<T> > result;

  if (sz == k)
  {
    std::list<T> l;

    for (int i = 0; i < sz; ++i) l.push_back(a[i]);

    result.push_back(l);

    return result;
  }
  else if (sz < k)
  {
    return result;
  }
  else if (k == 0)
  {
    std::list<T> l;
    result.push_back(l);

    return result;
  }

  for (int i = 0; i < sz; ++i)
  {
    T* nexta = (T*) malloc(sizeof(T) * (sz - 1));
    int ind = 0;

    for (int j = 0; j < sz; ++j)
    {
      if (j != i && a[j] > a[i]) nexta[ind++] = a[j];
    }

    std::list<std::list<T> > tmp_result = choose(nexta, ind, k - 1);
    typename std::list<std::list<T> >::iterator it = tmp_result.begin();
    for ( ; it != tmp_result.end(); ++it)
    {
      std::list<T> next = *it;
      next.push_back(a[i]);
      result.push_back(next);
    }

    free(nexta);
  }

  return result;
}

template <typename Graph>
class hash_real_edges {
public:
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef xmt_hash_table<int64_t, int64_t> hash_table_t;
  typedef xmt_hash_table<int64_t, double> dhash_table_t;

  hash_real_edges(Graph& gg, vertex_id_map<Graph>& vm, hash_table_t& tpc,
                  hash_table_t& eids, double* ons, double* gons,
                  dhash_table_t& el, dhash_table_t& er,
                  double* ew, int thresh) :
    g(gg), vid_map(vm), eid_map(get(_edge_id_map, gg)), real_edge_tp_count(tpc),
    real_edge_ids(eids), en_count(ons), good_en_count(gons), e_left(el),
    e_right(er), e_weight(ew), threshold(thresh), order(num_vertices(gg)) {}

  void operator()(edge_descriptor e)
  {
    vertex_descriptor v1 = source(e, g);
    vertex_descriptor v2 = target(e, g);
    int v1id = get(vid_map, v1);
    int v2id = get(vid_map, v2);
    int eid = get(eid_map, e);

    double wgt = e_weight[eid];
    en_count[eid] = 0;
    mt_incr(good_en_count[eid], wgt);

#ifdef DEBUG
    std::cout << "init good_en_count[" << eid << " (" << v1id << ", " << v2id
              << ")]: " << wgt << std::endl;
#endif

    order_pair(v1id, v2id);

    int64_t key = v1id * order + v2id;

#ifdef DEBUG
    std::cout << "hash_real_edges: " << eid << " (" << v1id << ", " << v2id
              << "): " << key << std::endl;
#endif

    real_edge_tp_count.insert(key, 0);
    real_edge_ids.insert(key, eid);

#ifdef RECTANGLES
    e_left.insert(key, 0);
    e_right.insert(key, 0);
#endif
  }

private:
  Graph& g;
  vertex_id_map<Graph>& vid_map;
  edge_id_map<Graph> eid_map;
  hash_table_t& real_edge_tp_count;
  hash_table_t& real_edge_ids;
  double* en_count;
  double* good_en_count;
  dhash_table_t& e_left;
  dhash_table_t& e_right;
  double* e_weight;
  int threshold;
  uint64_t order;
};

typedef int64_t intmax_t;

class hashvis {
public:
  hashvis() {}

  void operator()(int64_t key, int64_t value)
  {
    std::cout << "Table[" << key << "]: " << value << std::endl << std::flush;
  }
};

class dhashvis {
public:
  dhashvis() {}

  void operator()(int64_t key, double value)
  {
    std::cout << "Table[" << key << "]: " << value << std::endl << std::flush;
  }
};

template <typename Graph>
class fake_edge_finder {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
  typedef xmt_hash_table<int64_t, int64_t> hash_table_t;
  typedef xmt_hash_table<int64_t, double> dhash_table_t;

  fake_edge_finder(Graph& gg, hash_table_t& real, hash_table_t& fake,
                   dhash_table_t& el, dhash_table_t& er, int thresh = THRESH) :
    g(gg), threshold(thresh), next(0), real_edges(real),
    fake_edges(fake), e_left(el), e_right(er), eid_map(get(_edge_id_map, gg)) {}

  bool visit_test(vertex_descriptor src)
  {
    my_degree = out_degree(src, g);
    return (my_degree <= threshold);
  }

  void pre_visit(vertex_descriptor src)
  {
    my_edges = (size_type*) malloc(sizeof(size_type) * my_degree);
  }

  void operator()(edge_descriptor e, vertex_descriptor src,
                  vertex_descriptor dest)
  {
    int my_ind = mt_incr(next, 1);
    my_edges[my_ind] = get(eid_map, e);
  }

  void post_visit(vertex_descriptor src)
  {
    std::list<int> l;
    typedef std::list<std::list<size_type> > list_list_t;
    typedef typename std::list<size_type>::iterator l_iter_t;
    typedef typename std::list<std::list<size_type> >::iterator ll_iter_t;

    ll_iter_t it;
    vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);
    edge_iterator edgs = edges(g);

    unsigned int order = num_vertices(g);

    std::list<std::list<size_type> > all_pairs =
      choose(my_edges, static_cast<int>(my_degree), 2);

    it = all_pairs.begin();
    for ( ; it != all_pairs.end(); ++it)
    {
      std::list<size_type> l = *it;
      l_iter_t lit = l.begin();
      size_type eid1 = *lit++;
      size_type eid2 = *lit++;
      edge_descriptor e1 = edgs[eid1];
      edge_descriptor e2 = edgs[eid2];

      vertex_descriptor v1 = src == source(e1, g) ? target(e1, g) :
                                                    source(e1, g);
      vertex_descriptor v2 = src == source(e2, g) ? target(e2, g) :
                                                    source(e2, g);
      int64_t v1id = get(vid_map, v1);
      int64_t v2id = get(vid_map, v2);

      order_pair(v1id, v2id);

      int64_t key = v1id * order + v2id;

      if (!real_edges.member(key))
      {
        if (!fake_edges.member(key))
        {
          fake_edges.insert(key, 1);
          e_left.insert(key, 0);
          e_right.insert(key, 0);
        }
      }
    }

    free(my_edges);
  }

  Graph& g;
  int threshold;
  int64_t my_degree;
  int next;
  size_type* my_edges;
  hash_table_t& real_edges;
  hash_table_t& fake_edges;
  dhash_table_t& e_left;
  dhash_table_t& e_right;
  edge_id_map<Graph> eid_map;
};

//
// error_accumulator_full
//
// This was used to compute corrections when we were counting the
// edge 1-neighborhoods.  Ultimately, we didn't use this for the
// resolution limit paper.  If we revisit edge 1-neighborhoods,
// this deprecated code will have to be updated to accommodate the
// change in eweight from hashtable to array.  The pairs of neighbors
// will have to be edge id's rather than target vertices.  See
// fake_edge_finder above for a reference.
template <typename Graph>
class error_accumulator_full {
public:
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef xmt_hash_table<int64_t, int64_t> hash_table_t;
  typedef xmt_hash_table<int64_t, double> dhash_table_t;

  error_accumulator_full(Graph& gg, hash_table_t& real, hash_table_t& fake,
                         dhash_table_t& ew, dhash_table_t& el,
                         dhash_table_t& er, int thresh = THRESH) :
    g(gg), threshold(thresh), next(0), real_edges(real), fake_edges(fake),
    eweight(ew), order(num_vertices(gg)), e_left(el), e_right(er)
  {
    eid_map = get(_edge_id_map, g);
    vid_map = get(_vertex_id_map, g);
  }

  bool visit_test(vertex_descriptor src)
  {
    my_degree = out_degree(src, g);
    return (my_degree <= threshold);
  }

  void pre_visit(vertex_descriptor src)
  {
    my_id = get(vid_map, src);
    my_neighbors =
      (vertex_descriptor*) malloc(sizeof(vertex_descriptor) * my_degree);
  }

  void operator()(edge_descriptor e, vertex_descriptor src,
                  vertex_descriptor dest)
  {
    int64_t v1id = get(vid_map, src);
    int64_t v2id = get(vid_map, dest);
    int my_ind = mt_incr(next, 1);
    my_neighbors[my_ind] = dest;
    order_pair(v1id, v2id);
  }

  void post_visit(vertex_descriptor src)
  {
    typedef std::list<std::list<vertex_descriptor> > list_list_t;
    typedef typename std::list<vertex_descriptor>::iterator l_iter_t;
    typedef typename std::list<std::list<vertex_descriptor> >::iterator
            ll_iter_t;

    ll_iter_t it;
    unsigned int order = num_vertices(g);

    std::list<std::list<vertex_descriptor> > all_pairs =
      choose(my_neighbors, my_degree, 2);

    it = all_pairs.begin();
    for ( ; it != all_pairs.end(); ++it)
    {
      int vid = my_id;
      std::list<vertex_descriptor> l = *it;
      l_iter_t lit = l.begin();

      vertex_descriptor first = *lit++;
      vertex_descriptor second = *lit++;
      int v1id = get(vid_map, first);
      int v2id = get(vid_map, second);

      order_pair(v1id, v2id);

      int64_t key = v1id * order + v2id;

      order_pair(vid, v1id);

      int64_t l_key = vid * order + v1id;
      double wgt = eweight[l_key];

      e_left.update<hash_mt_incr<double> >(key, wgt, hash_mt_incr<double>());
//      e_left[key] += wgt;
      vid = my_id;

      order_pair(vid, v2id);

      int64_t r_key = vid * order + v2id;
      wgt = eweight[r_key];
      e_right.update<hash_mt_incr<double> >(key, wgt, hash_mt_incr<double>());
//      e_right[key] += wgt;
    }

    free(my_neighbors);
  }

  Graph& g;
  int threshold;
  int my_degree;
  int my_id;
  int next;
  unsigned int order;
  vertex_descriptor* my_neighbors;
  hash_table_t& real_edges;
  hash_table_t& fake_edges;
  dhash_table_t& e_left;
  dhash_table_t& e_right;
  dhash_table_t& eweight;
  vertex_id_map<Graph> vid_map;
  edge_id_map<Graph> eid_map;
};

// This version supports counting spokes only, so we don't really need the
// edge weights.  We'll just add 1 to each fake edge for each set of tent
// poles.  The distinction between "left" and "right" errors is deprecated,
// but the structures are retained to support future work on edge
// 1-neighborhoods.
template <typename Graph>
class error_accumulator {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
  typedef xmt_hash_table<int64_t, int64_t> hash_table_t;
  typedef xmt_hash_table<int64_t, double> dhash_table_t;

  error_accumulator(Graph& gg, hash_table_t& fake, double* ew,
                    dhash_table_t& el, dhash_table_t& er, int thresh = THRESH) :
    g(gg), threshold(thresh), next(0), fake_edges(fake),
    order(num_vertices(gg)), eweight(ew), e_left(el), e_right(er)
  {
    eid_map = get(_edge_id_map, g);
    vid_map = get(_vertex_id_map, g);
  }

  bool visit_test(vertex_descriptor src)
  {
    my_degree = out_degree(src, g);
    return (my_degree <= threshold);
  }

  void pre_visit(vertex_descriptor src)
  {
    my_id = get(vid_map, src);
    my_edges = (size_type*) malloc(sizeof(size_type) * my_degree);
  }

  void operator()(edge_descriptor e, vertex_descriptor src,
                  vertex_descriptor dest)
  {
    int my_ind = mt_incr(next, 1);
    my_edges[my_ind] = get(eid_map, e);
  }

  void post_visit(vertex_descriptor src)
  {
    typedef std::list<std::list<size_type> > list_list_t;
    typedef typename std::list<size_type>::iterator l_iter_t;
    typedef typename std::list<std::list<size_type> >::iterator ll_iter_t;

    edge_iterator edgs = edges(g);

    ll_iter_t it;
    unsigned int order = num_vertices(g);

    std::list<std::list<size_type> > all_pairs =
      choose(my_edges, my_degree, 2);

    it = all_pairs.begin();
    for ( ; it != all_pairs.end(); ++it)
    {
      std::list<size_type> l = *it;
      l_iter_t lit = l.begin();

      size_type eid1 = *lit++;
      size_type eid2 = *lit++;
      edge_descriptor e1 = edgs[eid1];
      edge_descriptor e2 = edgs[eid2];

      vertex_descriptor v1 = src == source(e1, g) ? target(e1, g) :
                                                    source(e1, g);
      vertex_descriptor v2 = src == source(e2, g) ? target(e2, g) :
                                                    source(e2, g);
      size_type v1id = get(vid_map, v1);
      size_type v2id = get(vid_map, v2);

      if (v2id < v1id) // The lower id is considered "left".
      {
        size_type tmp = eid1;
        eid1 = eid2;
        eid2 = tmp;
      }

      order_pair(v1id, v2id);
      int64_t key = v1id * order + v2id;
      if (!fake_edges.member(key)) continue;


      double l_wgt = eweight[eid1];
      double r_wgt = eweight[eid2];

      e_left.update<hash_mt_incr<double> >(key, l_wgt, hash_mt_incr<double>());
      e_right.update<hash_mt_incr<double> >(key, r_wgt, hash_mt_incr<double>());

      double lcur, rcur;
      e_left.lookup(key, lcur);
      e_right.lookup(key, rcur);
//      e_left[key] += 1;  // Simply note how many tent poles.
//      e_right[key] += 1;
    }

    free(my_edges);
  }

  Graph& g;
  int threshold;
  int next;
  hash_table_t& fake_edges;
  unsigned int order;
  double* eweight;
  dhash_table_t& e_left;
  dhash_table_t& e_right;
  vertex_id_map<Graph> vid_map;
  edge_id_map<Graph> eid_map;
  size_type* my_edges;
  int my_degree;
  int my_id;
};

template <typename Graph>
class ncount_correcter {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
  typedef xmt_hash_table<int64_t, int64_t> hash_table_t;
  typedef xmt_hash_table<int64_t, double> dhash_table_t;

  ncount_correcter(Graph& gg, hash_table_t& real, double* enc, double* genc,
                   double* ew, dhash_table_t& el, dhash_table_t& er,
                   int thresh = THRESH) :
    g(gg), order(num_vertices(gg)), vid_map(get(_vertex_id_map, g)),
    eid_map(get(_edge_id_map, g)), threshold(thresh), next(0),
    en_count(enc), good_en_count(genc), eweight(ew), real_edges(real),
    e_left(el), e_right(er) {}

  bool visit_test(vertex_descriptor src)
  {
    my_degree = out_degree(src, g);
    return (my_degree <= threshold);
  }

  void pre_visit(vertex_descriptor src)
  {
    my_id = get(vid_map, src);
    my_edges = (size_type*) malloc(sizeof(size_type) * my_degree);
  }

  void operator()(edge_descriptor e, vertex_descriptor src,
                  vertex_descriptor dest)
  {
    int my_ind = mt_incr(next, 1);
    my_edges[my_ind] = get(eid_map, e);
  }

  void post_visit(vertex_descriptor src)
  {
    typedef std::list<std::list<size_type> > list_list_t;
    typedef typename std::list<size_type>::iterator l_iter_t;
    typedef typename std::list<std::list<size_type> >::iterator ll_iter_t;

    edge_iterator edgs = edges(g);

    ll_iter_t it;

    std::list<std::list<size_type> > all_pairs =
      choose(my_edges, my_degree, 2);

    it = all_pairs.begin();
    for ( ; it != all_pairs.end(); ++it)
    {
      std::list<size_type> l = *it;
      l_iter_t lit = l.begin();

      size_type eid1 = *lit++;
      size_type eid2 = *lit++;
      edge_descriptor e1 = edgs[eid1];
      edge_descriptor e2 = edgs[eid2];

      vertex_descriptor v1 = src == source(e1, g) ? target(e1, g) :
                                                    source(e1, g);
      vertex_descriptor v2 = src == source(e2, g) ? target(e2, g) :
                                                    source(e2, g);
      size_type v1id = get(vid_map, v1);
      size_type v2id = get(vid_map, v2);
//      int v1id_bak = v1id;
//      int v2id_bak = v2id;

      if (v2id < v1id)   // The lower id is considered "left".
      {
        size_type tmp = eid1;
        eid1 = eid2;
        eid2 = tmp;
      }

      order_pair(v1id, v2id);

      uint64_t key = v1id * order + v2id;
      if (real_edges.member(key)) continue;

      double e1_wgt = eweight[eid1];
      double e2_wgt = eweight[eid2];
      double er = 0;
      double el = 0;

      if (real_edges.member(key))
      {
//        printf("%d -> %jd(%d,%d): add %lf to %jd, add %lf to %jd\n",
//               my_id, key, v1id_bak, v2id_bak, -0.25 * (er - r_wgt), l_key,
//               -0.25 * (el - l_wgt), r_key);
//        en_count[l_key] += (-0.25 * (er - r_wgt));
//        en_count[r_key] += (-0.25 * (el - l_wgt));
//        good_en_count[l_key] += (0.25 * (er - r_wgt));
//        good_en_count[r_key] += (0.25 * (el - l_wgt));
      }
      else
      {
        e_right.lookup(key, er);
        e_left.lookup(key, el);

//        printf("%d -> %jd(%d,%d): add %lf to %jd, add %lf to %jd\n",
//               my_id, key, v1id_bak, v2id_bak, 0.25 * (er - r_wgt), l_key,
//               0.25 * (el - l_wgt), r_key);
//        en_count[l_key] += (0.25 * (er - r_wgt));
//        en_count[r_key] += (0.25 * (el - l_wgt));
//        good_en_count[l_key] += (0.25 * (er - r_wgt));
//        good_en_count[r_key] += (0.25 * (el - l_wgt));

        if (el > e1_wgt || er > e2_wgt)  // >1 tent top
        {
//          printf("%d -> %jd(%d,%d): add %lf to %jd, add %lf to %jd\n",
//                 my_id, key, v1id_bak, v2id_bak, r_wgt, l_key, l_wgt, r_key);

          // SPOKES.
          good_en_count[eid1] += e2_wgt;
          good_en_count[eid2] += e1_wgt;

#ifdef DEBUG
          std::cout << "SPOKE*: good_en_count[" << eid1 << " (" << my_id << ", "
                    << v1id << ")] (for " << eid2 << " (" << my_id  << ", "
                    << v2id << "): " << good_en_count[eid1] << std::endl;
          std::cout << "SPOKE*: good_en_count[" << eid2 << " (" << my_id << ", "
                    << v2id << ")] (for " << eid2 << " (" << my_id  << ", "
                    << v1id << "): " << good_en_count[eid2] << std::endl;
#endif
        }
      }
    }

    free(my_edges);
  }

  Graph& g;
  size_type order;
  vertex_id_map<Graph> vid_map;
  edge_id_map<Graph> eid_map;
  int threshold;
  int next;
  double* en_count;
  double* good_en_count;
  double* eweight;
  hash_table_t& real_edges;
  dhash_table_t& e_left;
  dhash_table_t& e_right;
  int my_degree;
  int my_id;
  size_type* my_edges;
};

template <typename Graph>
class goodness_correcter {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef xmt_hash_table<int64_t, double> dhash_table_t;

  goodness_correcter(Graph& gg, double* vnc, dhash_table_t& enc,
                     dhash_table_t& genc, dhash_table_t& ew,
                     int thresh = THRESH) :
    g(gg), threshold(thresh), en_count(enc), vn_count(vnc),
    good_en_count(genc), eweight(ew), order(num_vertices(gg))
  {
    vid_map = get(_vertex_id_map, g);
  }

  void operator()(vertex_descriptor src, vertex_descriptor dest)
  {
    size_type my_degree = out_degree(src, g);

    if (my_degree > threshold) return;

    size_type dest_deg = out_degree(dest, g);

    if (dest_deg > threshold) return;

    size_type src_id = get(vid_map, src);
    size_type dest_id = get(vid_map, dest);

    if (src_id > dest_id) return;

    int64_t key = src_id * order + dest_id;
    double my_weight = eweight[key];
    double num_bad = vn_count[src_id] + vn_count[dest_id]
                     - my_weight - good_en_count[key];

    good_en_count[key] = en_count[key] - num_bad;
  }

  void post_visit() {}

private:
  Graph& g;
  int threshold;
  size_type order;
  double* vn_count;
  dhash_table_t& en_count;
  dhash_table_t& good_en_count;
  dhash_table_t& eweight;
  vertex_id_map<Graph> vid_map;
};

template <typename Graph>
class v_one_neighborhood_count : public default_triangles_visitor<Graph> {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef xmt_hash_table<int64_t, int64_t> hash_table_t;
  typedef xmt_hash_table<int64_t, double> dhash_table_t;

  v_one_neighborhood_count(Graph& gg, double* ct, vertex_id_map<Graph>& vm,
                           double* ew, double* enc, double* genc) :
    g(gg), order(num_vertices(gg)), vipm(vm), count(ct), eweight(ew),
    en_count(enc), good_en_count(genc) {}

/*
  void operator()(vertex_descriptor v1, vertex_descriptor v2,
                  vertex_descriptor v3)
  {
    size_type v1id = get(vipm, v1);
    size_type v2id = get(vipm, v2);
    size_type v3id = get(vipm, v3);

    // Can't have this when size_type is equal to vertex_descriptor.
    this->operator()(v1id, v2id, v3id);
  }
*/
  void operator()(size_type e1id, size_type e2id, size_type e3id)
  {
    mt_incr(good_en_count[e1id], eweight[e2id] + eweight[e3id]);
    mt_incr(good_en_count[e2id], eweight[e1id] + eweight[e3id]);
    mt_incr(good_en_count[e3id], eweight[e1id] + eweight[e2id]);

#ifdef DEBUG
    std::cout << "v_one_neigh::op(): good_en_count[" << e1id << "]: "
              << good_en_count[e1id] << std::endl;
    std::cout << "v_one_neigh::op(): good_en_count[" << e2id << "]: "
              << good_en_count[e2id] << std::endl;
    std::cout << "v_one_neigh::op(): good_en_count[" << e3id << "]: "
              << good_en_count[e3id] << std::endl;
#endif
  }

private:
  Graph& g;
  size_type order;
  vertex_id_map<Graph>& vipm;
  double* count;
  double* eweight;
  double* en_count;
  double* good_en_count;
};

/*! ********************************************************************
    \class neighbor_counts
    \brief
 *  ********************************************************************
*/
template <typename Graph>
class neighbor_counts {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
  typedef typename graph_traits<Graph>::out_edge_iterator out_edge_iterator;
  typedef xmt_hash_table<int64_t, int64_t> hash_table_t;
  typedef xmt_hash_table<int64_t, double> dhash_table_t;

  neighbor_counts(Graph& gg, double* ew, int thr = THRESH) :
    g(gg), eweight(ew), thresh(thr) {}

  void run(double* vn_count, double* en_count,
           double* good_en_count)
  {
    vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);
    edge_id_map<Graph> eid_map = get(_edge_id_map, g);

    unsigned int order = num_vertices(g);
    int size = num_edges(g);
    hashvis h;
    dhashvis dh;

#ifdef DEBUG
    if (order < 50) print(g);
#endif

    vertex_iterator verts = vertices(g);

    memset(vn_count, 0, order * sizeof(double));

#ifdef DEBUG
    double startTime = getWtime();
#endif

    #pragma mta assert nodep
    for (unsigned int i = 0; i < order; ++i)
    {
      vertex_descriptor v = verts[i];
      int deg = out_degree(v, g);
      out_edge_iterator eiter = out_edges(v, g);

      #pragma mta assert nodep
      for (int j = 0; j < deg; ++j)
      {
        edge_descriptor e = eiter[j];
        size_type eid = get(eid_map, e);
        mt_incr(vn_count[i], eweight[eid]);
      }
    }

#ifdef DEBUG
    double endTime = getWtime();
    std::cout << "Time: " << endTime - startTime << std::endl;
#endif

    hash_table_t real_edges(static_cast<hash_table_t::size_type>(1.5 * size));
    hash_table_t real_edge_ids(static_cast<hash_table_t::size_type>(1.5 *
                                                                    size));

#ifdef RECTANGLES
    // Need estimat. num.
    hash_table_t fake_edges(static_cast<hash_table_t::size_type>(FE_FACTOR *
                                                                 size));

    // Need to est. num.
    dhash_table_t e_left(static_cast<dhash_table_t::size_type>(FE_FACTOR *
                                                               size));

    // Fake edges.
    dhash_table_t e_right(static_cast<dhash_table_t::size_type>(FE_FACTOR *
                                                                size));

    hash_real_edges<Graph> hre(g, vid_map, real_edges, real_edge_ids,
                               en_count, good_en_count, e_left, e_right,
                               eweight, thresh);
#else
    dhash_table_t e_left(static_cast<dhash_table_t::size_type>(1));
    dhash_table_t e_right(static_cast<dhash_table_t::size_type>(1));

    hash_real_edges<Graph> hre(g, vid_map, real_edges, real_edge_ids,
                               en_count, good_en_count, e_left, e_right,
                               eweight, thresh);
#endif

    visit_edges(g, hre);

    v_one_neighborhood_count<Graph>
      onec(g, vn_count, vid_map, eweight, en_count, good_en_count);

    find_triangles(g, onec);

    edge_iterator edgs = edges(g);
    size_type m = num_edges(g);

    #pragma mta assert nodep
    for (size_type i = 0; i < m; ++i)
    {
      edge_descriptor e = edgs[i];
      size_type eid = get(eid_map, e);  // Should use property map.

      vertex_descriptor v1 = source(e, g);
      vertex_descriptor v2 = target(e, g);
      size_type v1id = get(vid_map, v1);
      size_type v2id = get(vid_map, v2);

      double vn1 = vn_count[v1id];
      double vn2 = vn_count[v2id];
      double incr = ((vn1 + vn2) - eweight[eid]);

      mt_incr(en_count[eid], incr);

      order_pair(v1id, v2id);
    }

#ifdef RECTANGLES
    fake_edge_finder<Graph> tpc(g, real_edges, fake_edges, e_left,
                                e_right, thresh);
    visit_adj(g, tpc);

#ifdef DEBUG
    std::cout << "real edges: " << real_edges.size() << std::endl 
              << "fake edges: " << fake_edges.size() << std::endl;
#endif

    error_accumulator<Graph> ea(g, fake_edges, eweight, e_left, e_right,
                                thresh);
    visit_adj(g, ea);

    ncount_correcter<Graph> ncc(g, real_edges, en_count, good_en_count,
                                eweight, e_left, e_right, thresh);
    visit_adj(g, ncc);

//    size_type* prefix_counts = 0;
//    size_type* started_nodes = 0;
//    size_type* num_threads;
//    memset(vn_count, 0, sizeof(double) * order);
//    goodness_correcter<Graph>
//      gc(g, vn_count, en_count, good_en_count, eweight, thresh);
//    visit_adj(g, gc, prefix_counts, started_nodes, num_threads);
//    printf("called goodness_corr\n");
//    fflush(stdout);
#endif
  }

private:
  Graph& g;
  double* eweight;
  int thresh;
};

template <typename Graph>
class weight_by_neighborhoods {
public:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
  typedef xmt_hash_table<int64_t, int64_t> hash_table_t;
  typedef xmt_hash_table<int64_t, double> dhash_table_t;

  weight_by_neighborhoods(Graph& gg, double* ew, int thr = THRESH) :
    g(gg), eweight(ew), thresh(thr) {}

  void run(int num_iter)
  {
    vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);
    edge_id_map<Graph> eid_map = get(_edge_id_map, g);
    unsigned int order = num_vertices(g);
    int size = num_edges(g);

    edge_iterator edgs = edges(g);

    for (int iter = 0; iter < num_iter; ++iter)
    {
      double* vn_count = (double*) calloc(order, sizeof(double));
      double* en_count = (double*) calloc(size, sizeof(double));
      double* good_en_count = (double*) calloc(size, sizeof(double));

      neighbor_counts<Graph> nc(g, eweight, thresh);
      nc.run(vn_count, en_count, good_en_count);

      for (int i = 0; i < size; ++i)
      {
        edge_descriptor e = edgs[i];
        size_type eid = get(eid_map, e);
        eweight[eid] = good_en_count[eid] / en_count[eid];

/*
        vertex_descriptor v1 = source(e, g);
        vertex_descriptor v2 = target(e, g);
        printf("eweight[%lu (%lu,%lu)]: (%f/%f) = %f\n", eid,
               get(vid_map, v1), get(vid_map, v2), good_en_count[eid],
               en_count[eid], good_en_count[eid] / en_count[eid]);
        fflush(stdout);
*/
      }

      free(vn_count);
      free(en_count);
      free(good_en_count);
    }
  }

private:
  Graph& g;
  double* eweight;
  int thresh;
};

}

#undef THRESH
#undef FE_FACTOR

#endif
