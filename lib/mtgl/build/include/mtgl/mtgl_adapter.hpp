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
/*! \file mtgl_adapter.hpp

    \brief This file contains common code used by adapters over graph
           data structures.

    \author Jon Berry (jberry@sandia.gov)

    \date 2007
*/
/****************************************************************************/

#ifndef MTGL_MTGL_ADAPTER_HPP
#define MTGL_MTGL_ADAPTER_HPP

#include <cstdio>
#include <iostream>
#include <fstream>

#include <mtgl/graph_traits.hpp>
#include <mtgl/mmap_traits.hpp>
#include <mtgl/mtgl_boost_property.hpp>
#include <mtgl/util.hpp>

// There are pragmas in this code specifically for the MTA
// The following line will disable the warning
#ifndef __MTA__
#pragma warning( disable: 4068 )
#endif

/// Namespace containing all the MTGL data structures and functions.
namespace mtgl {

enum vertex_id_map_t { _vertex_id_map = 0 };
enum edge_id_map_t { _edge_id_map = 1 };

template <typename Graph>
class vertex_id_map :
  public put_get_helper<int, vertex_id_map<Graph> > {
public:
  typedef typename graph_traits<Graph>::vertex_descriptor key_type;
  typedef typename graph_traits<Graph>::size_type value_type;
  vertex_id_map() {}

  value_type operator[] (const key_type& k) const { return k.get_id();  }
  value_type operator[] (const key_type* k) const { return k->get_id(); }
};

template <typename Graph>
class edge_id_map :
  public put_get_helper<int, edge_id_map<Graph> > {
public:
  typedef typename graph_traits<Graph>::edge_descriptor key_type;
  typedef typename graph_traits<Graph>::size_type value_type;
  edge_id_map() {}

  value_type operator[] (const key_type& k) const { return k.get_id();  }
  value_type operator[] (const key_type* k) const { return k->get_id(); }
};

template <typename Graph>
inline vertex_id_map<Graph> get(vertex_id_map_t, Graph& ga)
{
  return vertex_id_map<Graph>();
}

template <typename Graph>
inline edge_id_map<Graph> get(edge_id_map_t, Graph& ga)
{
  return edge_id_map<Graph>();
}

template <typename Graph, typename T>
class vertex_property_map :
  public array_property_map<T, vertex_id_map<Graph> > {
public:
  vertex_property_map(Graph& g) :
    array_property_map<T, vertex_id_map<Graph> >(
      (T*) malloc(num_vertices(g) * sizeof(T)), get(_vertex_id_map, g)) {}

  ~vertex_property_map() { free(this->data); }
};

template <typename Graph, typename T>
class edge_property_map :
  public array_property_map<T, edge_id_map<Graph> > {
public:
  edge_property_map(Graph& g) :
    array_property_map<T, edge_id_map<Graph> >(
      (T*) malloc(num_edges(g) * sizeof(T)), get(_edge_id_map, g)) {}

  ~edge_property_map() { free(this->data); }
};

struct keep_all {
  template <typename T>
  bool operator()(const T&) const { return true; }
};

template <typename Graph, typename visitor>
void
visit_edges(Graph& ga,
            typename graph_traits<Graph>::vertex_descriptor v,
            visitor& vis,
            typename graph_traits<Graph>::size_type par_cutoff = 5000,
            bool use_future = true,
            int directed = DIRECTED)
{
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::out_edge_iterator out_edge_iterator;
  typedef typename graph_traits<Graph>::size_type size_type;

  out_edge_iterator eit = out_edges(v, ga);
  size_type deg = out_degree(v, ga);

  if (deg >= par_cutoff)
  {
    if (use_future)
    {
      #pragma mta assert parallel
      #pragma mta loop future
      for (size_type i = 0; i < deg; i++)
      {
        edge_descriptor e = eit[i];
        vis(e);
      }
    }
    else
    {
      #pragma mta assert parallel
      for (size_type i = 0; i < deg; i++)
      {
        edge_descriptor e = eit[i];
        vis(e);
      }
    }
  }
  else
  {
    for (size_type i = 0; i < deg; i++)
    {
      edge_descriptor e = eit[i];
      vis(e);
    }
  }
}

template <typename Graph, typename visitor>
void visit_edges_filtered(Graph& ga,
                          typename graph_traits<Graph>::vertex_descriptor v,
                          visitor& vis, int par_cutoff = 5000,
                          bool use_future = true)
{
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::out_edge_iterator out_edge_iterator;

  out_edge_iterator eit = out_edges(v, ga);
  size_type deg = out_degree(v, ga);

  if (deg >= par_cutoff)
  {
    if (use_future)
    {
      #pragma mta assert parallel
      #pragma mta loop future
      for (size_type i = 0; i < deg; i++)
      {
        edge_descriptor e = eit[i];

        if (vis.visit_test(e)) vis(e);
      }
    }
    else
    {
      #pragma mta assert parallel
      for (size_type i = 0; i < deg; i++)
      {
        edge_descriptor e = eit[i];

        if (vis.visit_test(e)) vis(e);
      }
    }
  }
  else
  {
    for (size_type i = 0; i < deg; i++)
    {
      edge_descriptor e = eit[i];

      if (vis.visit_test(e)) vis(e);
    }
  }
}

template <typename Graph, typename visitor>
void visit_edges(Graph& ga,
                 visitor vis, bool use_future = true)
{
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
  typedef typename graph_traits<Graph>::size_type size_type;

  edge_iterator edgs = edges(ga);
  size_type size = num_edges(ga);

  if (use_future)
  {
    #pragma mta assert parallel
    #pragma mta loop future
    for (size_type i = 0; i < size; i++)
    {
      edge_descriptor e = edgs[i];
      vis(e);
    }
  }
  else
  {
    #pragma mta assert parallel
    for (size_type i = 0; i < size; i++)
    {
      edge_descriptor e = edgs[i];
      vis(e);
    }
  }
}

template <typename Graph, typename visitor>
void visit_edges_filtered(Graph& ga,
                          visitor vis, bool use_future = true)
{
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
  typedef typename graph_traits<Graph>::size_type size_type;

  edge_iterator edgs = edges(ga);
  size_type size = num_edges(ga);

  if (use_future)
  {
    #pragma mta assert parallel
    #pragma mta loop future
    for (size_type i = 0; i < size; i++)
    {
      edge_descriptor e = edgs[i];

      if (vis.visit_test(e)) vis(e);
    }
  }
  else
  {
    #pragma mta assert parallel
    for (size_type i = 0; i < size; i++)
    {
      edge_descriptor e = edgs[i];

      if (vis.visit_test(e)) vis(e);
    }
  }
}

template <typename Graph>
void print(Graph& ga)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::out_edge_iterator out_edge_iterator;

  vertex_id_map<Graph> vipm = get(_vertex_id_map, ga);
  edge_id_map<Graph> eipm = get(_edge_id_map, ga);

  vertex_iterator verts = vertices(ga);
  size_type order = num_vertices(ga);
  for (size_type i = 0; i < order; ++i)
  {
    vertex_descriptor v = verts[i];
    size_type vid = get(vipm, v);

    std::cout << vid << " : { ";

    out_edge_iterator inc_edgs = out_edges(v, ga);
    size_type out_deg = out_degree(v, ga);
    for (size_type j = 0; j < out_deg; ++j)
    {
      edge_descriptor e = inc_edgs[j];
      size_type eid = get(eipm, e);

      vertex_descriptor src = source(e, ga);
      vertex_descriptor trg = target(e, ga);
      size_type sid = get(vipm, src);
      size_type tid = get(vipm, trg);

      std::cout << "{" << eid << "}(" << sid << ", " << tid << ") ";
    }

    std::cout << "}" << std::endl;
  }
}

template <typename Graph>
void print_graph_viz(FILE* stream, Graph& ga)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::out_edge_iterator out_edge_iterator;

  fprintf(stream, "digraph \"G\" {\n");
  fprintf(stream, "graph[overlap=false splines=true size=\"7.5,7.5\"]\n");

  vertex_id_map<Graph> vipm = get(_vertex_id_map, ga);
  vertex_iterator verts = vertices(ga);
  size_type order = num_vertices(ga);

  for (size_type i = 0; i < order; ++i)
  {
    vertex_descriptor v = verts[i];
    int vid = get(vipm, v);

    out_edge_iterator inc_edgs = out_edges(v, ga);
    size_type out_deg = out_degree(v, ga);
    for (size_type j = 0; j < out_deg; ++j)
    {
      edge_descriptor e = inc_edgs[j];

      vertex_descriptor src = source(e, ga);
      vertex_descriptor trg = target(e, ga);

      int sid = get(vipm, src);
      int tid = get(vipm, trg);

      if (sid == vid) fprintf(stream, "%d -> %d;\n", sid, tid);
    }
  }

  fprintf(stream, "};\n");
}

template <typename Graph>
void print_edges(Graph& ga)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;

  vertex_id_map<Graph> vipm = get(_vertex_id_map, ga);
  edge_id_map<Graph> eipm = get(_edge_id_map, ga);

  edge_iterator edgs = edges(ga);
  size_type num_edgs = num_edges(ga);
  for (size_type i = 0; i < num_edgs; ++i)
  {
    edge_descriptor e = edgs[i];
    size_type eid = get(eipm, e);

    vertex_descriptor src = source(e, ga);
    vertex_descriptor trg = target(e, ga);

    size_type sid = get(vipm, src);
    size_type tid = get(vipm, trg);

    std::cout << "{" << eid << "}(" << sid << ", " << tid << ")" << std::endl;
  }
}

template <typename Graph>
void print_adj(Graph& ga)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::adjacency_iterator adjacency_iterator;

  vertex_id_map<Graph> vipm = get(_vertex_id_map, ga);

  vertex_iterator verts = vertices(ga);
  size_type num_verts = num_vertices(ga);
  for (size_type i = 0; i < num_verts; ++i)
  {
    vertex_descriptor src = verts[i];
    size_type sid = get(vipm, src);

    std::cout << sid << " : { ";

    adjacency_iterator adj_verts = adjacent_vertices(src, ga);
    size_type out_deg = out_degree(src, ga);
    for (size_type j = 0; j < out_deg; ++j)
    {
      vertex_descriptor trg = adj_verts[j];
      size_type tid = get(vipm, trg);

      std::cout << "(" << sid << ", " << tid << ") ";
    }

    std::cout << "}" << std::endl;
  }
}

template <typename Graph>
void separate_by_degree(Graph& g, int*& indicesOfBigs, int& num_big,
                        int*& indicesOfSmalls, int& num_small,
                        int deg_threshold = 5000)
{
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::size_type size_type;

  num_big   = 0;
  num_small = 0;
  int max_deg = 0;

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);

  vertex_iterator viter = vertices(g);
  size_type ord = num_vertices(g);

  #pragma mta assert parallel
  for (size_type i = 0; i < ord; ++i)
  {
    vertex_descriptor v = viter[i];

    size_type deg = out_degree(v, g);

    if (deg > max_deg) max_deg = deg;

    if (deg > deg_threshold)
    {
      mt_incr(num_big, 1);
    }
    else
    {
      mt_incr(num_small, 1);
    }
  }

  indicesOfBigs = (int*) malloc(sizeof(int) * num_big);
  indicesOfSmalls = (int*) malloc(sizeof(int) * num_small);

  num_big = num_small = 0;

  #pragma mta assert parallel
  for (size_type i = 0; i < ord; ++i)
  {
    vertex_descriptor v = viter[i];
    int id = get(vid_map, v);

    if (out_degree(v, g) > deg_threshold)
    {
      int mine = mt_incr(num_big, 1);
      indicesOfBigs[mine] = id;
    }
    else
    {
      int mine = mt_incr(num_small, 1);
      indicesOfSmalls[mine] = id;
    }
  }

  printf("sep by deg: %d, %d (max: %d)\n", num_big, num_small, max_deg);
}

template <typename Graph>
void separate_by_in_degree(Graph& g, int*& indicesOfBigs, int& num_big,
                           int*& indicesOfSmalls, int& num_small,
                           int deg_threshold = 5000)
{
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::size_type size_type;

  num_big   = 0;
  num_small = 0;
  int max_deg = 0;

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);

  vertex_iterator viter = vertices(g);
  size_type ord = num_vertices(g);

  #pragma mta assert parallel
  for (size_type i = 0; i < ord; ++i)
  {
    vertex_descriptor v = viter[i];

    size_type deg = in_degree(v, g);

    if (deg > max_deg) max_deg = deg;

    if (deg > deg_threshold)
    {
      mt_incr(num_big, 1);
    }
    else
    {
      mt_incr(num_small, 1);
    }
  }

  indicesOfBigs = (int*) malloc(sizeof(int) * num_big);
  indicesOfSmalls = (int*) malloc(sizeof(int) * num_small);
  num_big = num_small = 0;

  #pragma mta assert parallel
  for (size_type i = 0; i < ord; ++i)
  {
    vertex_descriptor v = viter[i];
    int id = get(vid_map, v);

    if (in_degree(v, g) > deg_threshold)
    {
      int mine = mt_incr(num_big, 1);
      indicesOfBigs[mine] = id;
    }
    else
    {
      int mine = mt_incr(num_small, 1);
      indicesOfSmalls[mine] = id;
    }
  }

  printf("sep by in: %d, %d (max: %d)\n", num_big, num_small, max_deg);
}

//
// tython_viz_graph(g, fn)
//
//     Produce output convenient for Titan's "tython" viz interface
//
template <typename Graph>
void
tython_viz_graph(Graph& g, char *name)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);
  size_type size = num_edges(g);
  size_type order = num_vertices(g);

  char gfname[256];
  char vfname[256];
  sprintf(gfname, "%s.csv", name);
  sprintf(vfname, "%s.vertices.csv", name);

  std::ofstream gf(gfname);
  std::ofstream vf(vfname);
  vf << "id property" << std::endl;
  gf << "source target" << std::endl;

  for (size_type i = 0; i < order; ++i) vf << i << " 0" << std::endl;

  edge_iterator edgs = edges(g);

  for (size_type i = 0; i < size; ++i)
  {
    edge_descriptor e = edgs[i];

    vertex_descriptor u = source(e, g);
    vertex_descriptor v = target(e, g);
    size_type uid = get(vid_map, u);
    size_type vid = get(vid_map, v);

    printf("edge: %lu, %lu\n", uid, vid);

    gf << uid << " " << vid << std::endl;
  }

  gf.close();
  vf.close();
}

}

#endif
