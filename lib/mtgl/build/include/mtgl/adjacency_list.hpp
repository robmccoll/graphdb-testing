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
/*! \file adjacency_list.hpp

    \brief An adjacency list graph implementation.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 3/17/2008
*/
/****************************************************************************/

#ifndef MTGL_ADJACENCY_LIST_HPP
#define MTGL_ADJACENCY_LIST_HPP

#include <limits>
#include <iostream>
#include <string>

#include <mtgl/mtgl_adapter.hpp>
#include <mtgl/dynamic_array.hpp>
#include <mtgl/xmt_hash_set.hpp>

namespace mtgl {

namespace detail {

template <typename Graph>
class al_edge_adapter {
private:
  typedef typename Graph::Vertex Vertex;
  typedef typename Graph::Edge Edge;

public:
  al_edge_adapter() : first(NULL), second(NULL),
                      id((std::numeric_limits<unsigned long>::max)()) {}
  al_edge_adapter(Edge& e) : first(e.from), second(e.to), id(e.id) {}
  al_edge_adapter(Edge* e) : first(e->from), second(e->to), id(e->id) {}

  al_edge_adapter(Vertex* v1, Vertex* v2, unsigned long eid) :
    first(v1), second(v2), id(eid) {}

  bool operator==(const al_edge_adapter& rhs) const
  { return id == rhs.id; }

public:
  Vertex* first;
  Vertex* second;
  unsigned long id;
};

/***/

template <typename Graph>
class al_adjacency_iterator {
private:
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename Graph::Vertex Vertex;
  typedef typename Graph::Edge Edge;

public:
  al_adjacency_iterator() : v(0) {}
  al_adjacency_iterator(Vertex* vert) : v(vert) {}

  vertex_descriptor operator[](unsigned long p) const
  {
    Edge* e = v->edge_list[p];
    return v == e->from ? e->to : e->from;
  }

private:
  Vertex* v;
};

/***/

template <typename Graph>
class al_out_edge_iterator {
private:
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename Graph::Vertex Vertex;
  typedef typename Graph::Edge Edge;

public:
  al_out_edge_iterator() : v(0) {}
  al_out_edge_iterator(Vertex* vert) : v(vert) {}

  edge_descriptor operator[](unsigned long p) const
  {
    Edge* e = v->edge_list[p];
    Vertex* v2 = v == e->from ? e->to : e->from;
    return edge_descriptor(v, v2, e->id);
  }

private:
  Vertex* v;
};

/***/

template <typename Graph>
class al_thread_vertex_iterator {
private:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename Graph::Vertex Vertex;

public:
  al_thread_vertex_iterator() : vertex_list(0) {}

  al_thread_vertex_iterator(Vertex** vl, size_type p) : vertex_list(vl + p) {}

  al_thread_vertex_iterator(const al_thread_vertex_iterator& rhs) :
    vertex_list(rhs.vertex_list) {}

  al_thread_vertex_iterator& operator=(const al_thread_vertex_iterator& rhs)
  {
    vertex_list = rhs.vertex_list;
    return *this;
  }

  al_thread_vertex_iterator& operator++()
  {
    ++vertex_list;
    return *this;
  }

  al_thread_vertex_iterator& operator++(int)
  {
    al_thread_vertex_iterator temp(*this);

    ++vertex_list;
    return temp;
  }

  vertex_descriptor operator*() const { return *vertex_list; }

  bool operator==(const al_thread_vertex_iterator& rhs) const
  { return vertex_list == rhs.vertex_list; }
  bool operator!=(const al_thread_vertex_iterator& rhs) const
  { return vertex_list != rhs.vertex_list; }

private:
  Vertex** vertex_list;
};

/***/

template <typename Graph>
class al_thread_edge_iterator {
private:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename Graph::Edge Edge;

public:
  al_thread_edge_iterator() : edge_list(0) {}

  al_thread_edge_iterator(Edge** el, size_type p) : edge_list(el + p) {}

  al_thread_edge_iterator(const al_thread_edge_iterator& rhs) :
    edge_list(rhs.edge_list) {}

  al_thread_edge_iterator& operator=(const al_thread_edge_iterator& rhs)
  {
    edge_list = rhs.edge_list;
    return *this;
  }

  al_thread_edge_iterator& operator++()
  {
    ++edge_list;
    return *this;
  }

  al_thread_edge_iterator& operator++(int)
  {
    al_thread_edge_iterator temp(*this);

    ++edge_list;
    return temp;
  }

  edge_descriptor operator*() const { return *edge_list; }

  bool operator==(const al_thread_edge_iterator& rhs) const
  { return edge_list == rhs.edge_list; }
  bool operator!=(const al_thread_edge_iterator& rhs) const
  { return edge_list != rhs.edge_list; }

private:
  Edge** edge_list;
};

/***/

template <typename Graph>
class al_thread_adjacency_iterator {
private:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename Graph::Vertex Vertex;
  typedef typename Graph::Edge Edge;

public:
  al_thread_adjacency_iterator() : v(0), pos(0) {}

  al_thread_adjacency_iterator(Vertex* vert, size_type p) : v(vert), pos(p) {}

  al_thread_adjacency_iterator(const al_thread_adjacency_iterator& rhs) :
    v(rhs.v), pos(rhs.pos) {}

  al_thread_adjacency_iterator&
  operator=(const al_thread_adjacency_iterator& rhs)
  {
    v = rhs.v;
    pos = rhs.pos;
    return *this;
  }

  al_thread_adjacency_iterator& operator++()
  {
    ++pos;
    return *this;
  }

  al_thread_adjacency_iterator& operator++(int)
  {
    al_thread_adjacency_iterator temp(*this);

    ++pos;
    return temp;
  }

  vertex_descriptor operator*() const
  {
    Edge* e = v->edge_list[pos];
    return v == e->from ? e->to : e->from;
  }

  bool operator==(const al_thread_adjacency_iterator& rhs) const
  { return v == rhs.v; }
  bool operator!=(const al_thread_adjacency_iterator& rhs) const
  { return v != rhs.v; }

private:
  Vertex* v;
  size_type pos;
};

/***/

template <typename Graph>
class al_thread_out_edge_iterator {
private:
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename Graph::Vertex Vertex;
  typedef typename Graph::Edge Edge;

public:
  al_thread_out_edge_iterator() : v(0), pos(0) {}

  al_thread_out_edge_iterator(Vertex* vert, size_type p) : v(vert), pos(p) {}

  al_thread_out_edge_iterator(const al_thread_out_edge_iterator& rhs) :
    v(rhs.v), pos(rhs.pos) {}

  al_thread_out_edge_iterator&
  operator=(const al_thread_out_edge_iterator& rhs)
  {
    v = rhs.v;
    pos = rhs.pos;
    return *this;
  }

  al_thread_out_edge_iterator& operator++()
  {
    ++pos;
    return *this;
  }

  al_thread_out_edge_iterator& operator++(int)
  {
    al_thread_out_edge_iterator temp(*this);

    ++pos;
    return temp;
  }

  edge_descriptor operator*() const
  {
    Edge* e = v->edge_list[pos];
    Vertex* v2 = v == e->from ? e->to : e->from;
    return edge_descriptor(v, v2, e->id);
  }

  bool operator==(const al_thread_out_edge_iterator& rhs) const
  { return v == rhs.v; }
  bool operator!=(const al_thread_out_edge_iterator& rhs) const
  { return v != rhs.v; }

private:
  Vertex* v;
  size_type pos;
};

}

/***/

template <typename DIRECTION = directedS>
class adjacency_list {
public:
  struct Vertex;
  struct Edge;

  typedef unsigned long size_type;
  typedef Vertex* vertex_descriptor;
  typedef detail::al_edge_adapter<adjacency_list> edge_descriptor;
  typedef Vertex** vertex_iterator;
  typedef detail::al_adjacency_iterator<adjacency_list> adjacency_iterator;
  typedef void in_adjacency_iterator;
  typedef Edge** edge_iterator;
  typedef detail::al_out_edge_iterator<adjacency_list> out_edge_iterator;
  typedef void in_edge_iterator;
  typedef detail::al_thread_vertex_iterator<adjacency_list>
          thread_vertex_iterator;
  typedef detail::al_thread_adjacency_iterator<adjacency_list>
          thread_adjacency_iterator;
  typedef void thread_in_adjacency_iterator;
  typedef detail::al_thread_edge_iterator<adjacency_list> thread_edge_iterator;
  typedef detail::al_thread_out_edge_iterator<adjacency_list>
          thread_out_edge_iterator;
  typedef void thread_in_edge_iterator;
  typedef DIRECTION directed_category;
  typedef vector_thread_iterators iterator_category;

  adjacency_list() : nVertices(0), nEdges(0) {}
  adjacency_list(const adjacency_list& g) { deep_copy(g); }
  ~adjacency_list() { clear(); }

  adjacency_list& operator=(const adjacency_list& rhs)
  {
    clear();
    deep_copy(rhs);
    return *this;
  }

  inline void clear();

  void init(size_type num_verts, size_type num_edges,
            size_type* sources, size_type* targets)
  {
    clear();
    addVertices(num_verts);
    addEdges(num_edges, sources, targets);
  }

  size_type get_order() const { return nVertices; }
  size_type get_size() const { return nEdges;  }

  size_type get_degree(vertex_descriptor v) const
  {
    return v->edge_list.size();
  }

  size_type get_out_degree(vertex_descriptor v) const
  {
    return v->edge_list.size();
  }

  vertex_iterator vertices() const { return vertex_list.get_data(); }

  adjacency_iterator adjacent_vertices(const vertex_descriptor& v) const
  {
    return adjacency_iterator(v);
  }

  edge_iterator edges() const { return edge_list.get_data(); }

  out_edge_iterator out_edges(const vertex_descriptor& v) const
  {
    return out_edge_iterator(v);
  }

  thread_vertex_iterator
  thread_vertices(size_type pos) const
  {
    return thread_vertex_iterator(vertex_list.get_data(), pos);
  }

  thread_adjacency_iterator
  thread_adjacent_vertices(const vertex_descriptor& v, size_type pos) const
  {
    return thread_adjacency_iterator(v, pos);
  }

  thread_edge_iterator
  thread_edges(size_type pos) const
  {
    return thread_edge_iterator(edge_list.get_data(), pos);
  }

  thread_out_edge_iterator
  thread_out_edges(const vertex_descriptor& v, size_type pos) const
  {
    return thread_out_edge_iterator(v, pos);
  }

  void print() const;
  void printAdj() const;
  void printDotfile(const std::string& filename) const;

  inline Vertex* addVertex();
  inline void addVertices(size_type num_verts);
  inline void removeVertices(size_type num_verts, Vertex** v);

  template <typename T>
  inline void removeVertices(size_type num_verts, size_type* v, 
                             edge_property_map<adjacency_list, T>* emap = 0);

  inline Edge* addEdge(Vertex* f, Vertex* t);
  Edge* addEdge(size_type f_id, size_type t_id)
  { return addEdge(vertex_list[f_id], vertex_list[t_id]); }

  inline void addEdges(size_type num_edges, Vertex** f, Vertex** t);
  inline void addEdges(size_type num_edges, size_type* f, size_type* t);
  inline void removeEdges(size_type num_edges, size_type* e);

  template <typename T>
  inline void removeEdges(size_type num_edges, Edge** e,
                          edge_property_map<adjacency_list, T>* emap = 0);

private:
  inline void deep_copy(const adjacency_list& rhs);

public:
  size_type nVertices;
  size_type nEdges;

  dynamic_array<Vertex*> vertex_list;
  dynamic_array<Edge*> edge_list;
};

/***/

template <>
class adjacency_list<bidirectionalS> {
public:
  adjacency_list()
  {
    fprintf(stderr, "** Error: Bidirectional direction not yet supported "
            "for adjacency_list. **\n");
    exit(1);
  }
};

/***/

template <typename DIRECTION>
struct adjacency_list<DIRECTION>::Edge {
  Edge(Vertex* f, Vertex* t, size_type i) : id(i), from(f), to(t) {}

  size_type id;
  Vertex* from;
  Vertex* to;
};

/***/

template <typename DIRECTION>
struct adjacency_list<DIRECTION>::Vertex {
  Vertex(size_type i) : id(i), num_edges_to_add(0) {}

  void removeEdge(Edge* e)
  {
    size_type lock_val = mt_readfe(num_edges_to_add);

    // Find the edge.
    size_type i = 0;
    for ( ; i < edge_list.size() && edge_list[i] != e; ++i);

    // If it isn't the last edge, replace the empty entry with the last edge.
    // This test catches the case when there is only one edge as, for this
    // case, edge_list.size() - 1 is 0.
    if (i < edge_list.size() - 1)
    {
      edge_list[i] = edge_list[edge_list.size() - 1];
    }

    // If the edge was found it is now in the last position.  Make the array
    // one smaller.
    if (i < edge_list.size()) edge_list.resize(edge_list.size() - 1);

    mt_write(num_edges_to_add, lock_val);
  }

  size_type id;
  size_type num_edges_to_add;
  dynamic_array<Edge*> edge_list;
};

/***/

template <typename DIRECTION>
void
adjacency_list<DIRECTION>::clear()
{
  nVertices = 0;
  nEdges = 0;

  #pragma mta assert parallel
  for (size_type i = 0; i < vertex_list.size(); ++i) delete vertex_list[i];

  #pragma mta assert parallel
  for (size_type i = 0; i < edge_list.size(); ++i) delete edge_list[i];

  // Empty the vertex and edge arrays.
  vertex_list.clear();
  edge_list.clear();
}

/***/

template <typename DIRECTION>
void
adjacency_list<DIRECTION>::print() const
{
  std::cout << std::endl << "Vertices: " << nVertices << std::endl;
  size_type nVerts = vertex_list.size();
  for (size_type i = 0; i < nVerts; ++i)
  {
    Vertex* v = vertex_list[i];

    std::cout << v->id << std::endl;
  }

  std::cout << std::endl << "Edges: " << nEdges << std::endl;
  size_type num_edges = edge_list.size();
  for (size_type i = 0; i < num_edges; ++i)
  {
    Edge* e = edge_list[i];

    std::cout << e->id << " (" << e->from->id << ", "
              << e->to->id << ")" << std::endl;
  }
}

/***/

template <typename DIRECTION>
void
adjacency_list<DIRECTION>::printAdj() const
{
  std::cout << std::endl << "Vertices: " << nVertices << std::endl;
  size_type nVerts = vertex_list.size();
  for (size_type i = 0; i < nVerts; ++i)
  {
    Vertex* v = vertex_list[i];

    std::cout << v->id << ":" << std::flush;

    for (size_type j = 0; j < v->edge_list.size(); ++j)
    {
      Edge* e = v->edge_list[j];

      std::cout << " " << e->id << " (" << e->from->id << ", "
                << e->to->id << ")" << std::flush;
    }

    std::cout << std::endl;
  }

  std::cout << std::endl;
}

/***/

template <typename DIRECTION>
void
adjacency_list<DIRECTION>::printDotfile(const std::string& filename) const
{
  FILE* GRAPHFILE = fopen(filename.c_str(), "w");
  fprintf(GRAPHFILE, "digraph gname {\n");
  fprintf(GRAPHFILE, "  ratio=.5;\n");
  fprintf(GRAPHFILE, "  margin=.1;\n\n");

  for (size_type i = 0; i < nEdges; ++i)
  {
    fprintf(GRAPHFILE, "  %7d -> %7d;\n", edge_list[i]->from->id,
            edge_list[i]->to->id);
  }

  fprintf(GRAPHFILE, "}\n");
  fclose(GRAPHFILE);
}

/***/

template <typename DIRECTION>
typename adjacency_list<DIRECTION>::Vertex*
adjacency_list<DIRECTION>::addVertex()
{
  size_type id = mt_readfe(nVertices);
  Vertex* v = new Vertex(id);
  vertex_list.push_back(v);
  mt_write(nVertices, nVertices + 1);

  return v;
}

/***/

template <typename DIRECTION>
void
adjacency_list<DIRECTION>::addVertices(size_type num_verts)
{
  size_type id = mt_readfe(nVertices);

  vertex_list.resize(id + num_verts);

  // I might be able to move the mt_write to here, but for now we need to
  // be conservative with unlocking.

  #pragma mta assert nodep
  for (size_type i = 0; i < num_verts; ++i)
  {
    vertex_list[id + i] = new Vertex(id + i);
  }

  mt_write(nVertices, id + num_verts);
}

/***/

template <typename DIRECTION>
void
adjacency_list<DIRECTION>::removeVertices(size_type num_verts, Vertex** v)
{
  size_type id = mt_readfe(nVertices);

  // Create a set that has all the vertices to be deleted.
  xmt_hash_set<size_type> delete_verts(2 * num_verts);

  #pragma mta assert parallel
  for (size_type i = 0; i < num_verts; ++i) delete_verts.insert(v[i]->id);

  // Count the number of edges in the graph that have one or more endpoints
  // that will be deleted.
  size_type total_edges = 0;
  for (size_type i = 0; i < id; ++i)
  {
    for (size_type j = 0; j < vertex_list[i]->edge_list.size(); ++j)
    {
      total_edges +=
        (delete_verts.member(vertex_list[i]->edge_list[j]->from->id) ||
         delete_verts.member(vertex_list[i]->edge_list[j]->to->id));
    }
  }

  // Loop over all the vertices and all their adjacency lists.  Mark an edge
  // to be deleted if either its source or target belong to the set of
  // vertices to be deleted.
  Edge** edges_to_delete = (Edge**) malloc(sizeof(Edge*) * total_edges);
  int edge_pos = 0;
  #pragma mta assert parallel
  for (size_type i = 0; i < id; ++i)
  {
    #pragma mta assert parallel
    for (size_type j = 0; j < vertex_list[i]->edge_list.size(); ++j)
    {
      if (!DIRECTION::is_directed() && 
          (vertex_list[i]->edge_list[j]->to->id == vertex_list[i]->id)) {
           continue;
      }
      if (delete_verts.member(vertex_list[i]->edge_list[j]->from->id) ||
          delete_verts.member(vertex_list[i]->edge_list[j]->to->id))
      {
        int cur_pos = mt_incr(edge_pos, 1);
        edges_to_delete[cur_pos] = vertex_list[i]->edge_list[j];
      }
    }
  }

  // Remove the edges and free the temporary memory.
  removeEdges<bool>(total_edges, edges_to_delete);
  free(edges_to_delete);

  // Remove the vertices from the vertex list that are in the area that is to
  // be deleted.
  #pragma mta assert nodep
  for (size_type i = 0; i < num_verts; ++i)
  {
    if (v[i]->id >= id - num_verts) vertex_list[v[i]->id] = 0;
  }

  // Remove the vertices from the vertex list that are not in the area that is
  // to be deleted.
  size_type moved_verts = 0;
  #pragma mta assert nodep
  for (size_type i = 0; i < num_verts; ++i)
  {
    if (v[i]->id < id - num_verts)
    {
      // Find the next last vertex in the area to be deleted that hasn't
      // already been deleted.
      size_type nmv = mt_incr(moved_verts, 1);
      while (vertex_list[id - 1 - nmv] == 0) nmv = mt_incr(moved_verts, 1);

      // Move the last vertex into the deleted vertex's position.
      vertex_list[v[i]->id] = vertex_list[id - 1 - nmv];

      // Change the moved vertex to have the id of the vertex it replaced.
      vertex_list[v[i]->id]->id = v[i]->id;
    }
  }

  vertex_list.resize(id - num_verts);

  mt_write(nVertices, id - num_verts);

  // Delete the memory associated with the deleted vertices.
  #pragma mta assert parallel
  for (size_type i = 0; i < num_verts; ++i) delete v[i];
}

/***/

template <typename DIRECTION>
template <typename T>
void
adjacency_list<DIRECTION>::removeVertices(
    size_type num_verts, size_type* v,
    edge_property_map<adjacency_list<DIRECTION>, T>* emap)
{
  size_type id = mt_readfe(nVertices);

  // Create a set that has all the vertices to be deleted.
  xmt_hash_set<size_type> delete_verts(2 * num_verts);

  #pragma mta assert parallel
  for (size_type i = 0; i < num_verts; ++i) delete_verts.insert(v[i]);

  // Count the number of edges in the graph that have one or more endpoints
  // that will be deleted.
  size_type total_edges = 0;
  for (size_type i = 0; i < id; ++i)
  {
    for (size_type j = 0; j < vertex_list[i]->edge_list.size(); ++j)
    {
      total_edges +=
        (delete_verts.member(vertex_list[i]->edge_list[j]->from->id) ||
         delete_verts.member(vertex_list[i]->edge_list[j]->to->id));
    }
  }

  // Loop over all the vertices and all their adjacency lists.  Mark an edge
  // to be deleted if either its source or target belong to the set of
  // vertices to be deleted.
  Edge** edges_to_delete = (Edge**) malloc(sizeof(Edge*) * total_edges);
  int edge_pos = 0;
  #pragma mta assert parallel
  for (size_type i = 0; i < id; ++i)
  {
    #pragma mta assert parallel
    for (size_type j = 0; j < vertex_list[i]->edge_list.size(); ++j)
    {
      if (!DIRECTION::is_directed() && 
          (vertex_list[i]->edge_list[j]->to->id == vertex_list[i]->id)) {
           continue;
      }
      if (delete_verts.member(vertex_list[i]->edge_list[j]->from->id) ||
          delete_verts.member(vertex_list[i]->edge_list[j]->to->id))
      {
        int cur_pos = mt_incr(edge_pos, 1);
        edges_to_delete[cur_pos] = vertex_list[i]->edge_list[j];
      }
    }
  }
  total_edges = edge_pos;

  // Remove the edges and free the temporary memory.
  removeEdges(total_edges, edges_to_delete, emap);
  free(edges_to_delete);

  // Remove the vertices from the vertex list that are in the area that is to
  // be deleted and delete their memory.
  #pragma mta assert nodep
  for (size_type i = 0; i < num_verts; ++i)
  {
    if (v[i] >= id - num_verts)
    {
      delete vertex_list[v[i]];
      vertex_list[v[i]] = 0;
    }
  }

  // Remove the vertices from the vertex list that are not in the area that is
  // to be deleted and delete their memory.
  size_type moved_verts = 0;
  #pragma mta assert nodep
  for (size_type i = 0; i < num_verts; ++i)
  {
    if (v[i] < id - num_verts)
    {
      // Find the next last vertex in the area to be deleted that hasn't
      // already been deleted.
      size_type nmv = mt_incr(moved_verts, 1);
      while (vertex_list[id - 1 - nmv] == 0) nmv = mt_incr(moved_verts, 1);

      delete vertex_list[v[i]];

      // Move the last vertex into the deleted vertex's position.
      vertex_list[v[i]] = vertex_list[id - 1 - nmv];

      // Change the moved vertex to have the id of the vertex it replaced.
      vertex_list[v[i]]->id = v[i];
    }
  }

  vertex_list.resize(id - num_verts);

  mt_write(nVertices, id - num_verts);
}

/***/

template <typename DIRECTION>
typename adjacency_list<DIRECTION>::Edge*
adjacency_list<DIRECTION>::addEdge(Vertex* f, Vertex* t)
{
  size_type id = mt_readfe(nEdges);
  Edge* e = new Edge(f, t, id);
  edge_list.push_back(e);

  f->edge_list.push_back(e);
  if (!DIRECTION::is_directed()) t->edge_list.push_back(e);

  mt_write(nEdges, nEdges + 1);

  return e;
}

/***/

template <typename DIRECTION>
void
adjacency_list<DIRECTION>::addEdges(size_type num_edges, Vertex** f, Vertex** t)
{
  size_type id = mt_readfe(nEdges);

  edge_list.resize(id + num_edges);

  // Add the new edges to the edges array.
  #pragma mta assert nodep
  for (size_type i = 0; i < num_edges; ++i)
  {
    edge_list[id + i] = new Edge(f[i], t[i], id + i);
  }

  // Count the number of edges to add to each vertex's edge list.
  #pragma mta assert nodep
  for (size_type i = 0; i < num_edges; ++i)
  {
    mt_incr(f[i]->num_edges_to_add, 1);
    if (!DIRECTION::is_directed()) mt_incr(t[i]->num_edges_to_add, 1);
  }

  // TODO: Instead of using a hash set and looping over the edges, I could
  //       also just loop over all the vertices in the graph.  Vertices that
  //       don't have edges being added will have num_edges_to_add set to
  //       0, so the reserve call won't do anything.

  // Resize each vertex's edge list and reset each vertex's num_edges_to_add
  // to 0.
  xmt_hash_set<size_type> unique_verts(2 * num_edges);

  #pragma mta assert parallel
  for (size_type i = 0; i < num_edges; ++i)
  {
    if (unique_verts.insert(f[i]->id).second)
    {
      f[i]->edge_list.reserve(f[i]->edge_list.size() + f[i]->num_edges_to_add);
      f[i]->num_edges_to_add = 0;
    }

    if (!DIRECTION::is_directed())
    {
      if (unique_verts.insert(t[i]->id).second)
      {
        t[i]->edge_list.reserve(t[i]->edge_list.size() +
                                t[i]->num_edges_to_add);
        t[i]->num_edges_to_add = 0;
      }
    }
  }

  // Add the new edges to each vertex's edge list.  This can be done
  // "unsafely" since we have already resized all the edge lists
  // appropriately.
  #pragma mta assert nodep
  for (size_type i = 0; i < num_edges; ++i)
  {
    f[i]->edge_list.unsafe_push_back(edge_list[id + i]);
    if (!DIRECTION::is_directed())
    {
      t[i]->edge_list.unsafe_push_back(edge_list[id + i]);
    }
  }

  mt_write(nEdges, id + num_edges);
}

/***/

template <typename DIRECTION>
void
adjacency_list<DIRECTION>::addEdges(size_type num_edges, size_type* f,
                                    size_type* t)
{
  size_type id = mt_readfe(nEdges);

  edge_list.resize(id + num_edges);

  // Add the new edges to the edges array.
  #pragma mta assert nodep
  for (size_type i = 0; i < num_edges; ++i)
  {
    edge_list[id + i] = new Edge(vertex_list[f[i]], vertex_list[t[i]], id + i);
  }

  // Count the number of edges to add to each vertex's edge list.
  #pragma mta assert nodep
  for (size_type i = 0; i < num_edges; ++i)
  {
    mt_incr(vertex_list[f[i]]->num_edges_to_add, 1);
    if (!DIRECTION::is_directed())
    {
      mt_incr(vertex_list[t[i]]->num_edges_to_add, 1);
    }
  }

  // TODO: Instead of using a hash set and looping over the edges, I could
  //       also just loop over all the vertices in the graph.  Vertices that
  //       don't have edges being added will have num_edges_to_add set to
  //       0, so the reserve call won't do anything.

  // Resize each vertex's edge list and reset each vertex's num_edges_to_add
  // to 0.
  xmt_hash_set<size_type> unique_verts(2 * num_edges);

  #pragma mta assert parallel
  for (size_type i = 0; i < num_edges; ++i)
  {
    if (unique_verts.insert(vertex_list[f[i]]->id).second)
    {
      vertex_list[f[i]]->edge_list.reserve(vertex_list[f[i]]->edge_list.size() +
                                           vertex_list[f[i]]->num_edges_to_add);
      vertex_list[f[i]]->num_edges_to_add = 0;
    }

    if (!DIRECTION::is_directed())
    {
      if (unique_verts.insert(vertex_list[t[i]]->id).second)
      {
        vertex_list[t[i]]->edge_list.reserve(
          vertex_list[t[i]]->edge_list.size() +
          vertex_list[t[i]]->num_edges_to_add);
        vertex_list[t[i]]->num_edges_to_add = 0;
      }
    }
  }

  // Add the new edges to each vertex's edge list.  This can be done
  // "unsafely" since we have already resized all the edge lists
  // appropriately.
  #pragma mta assert nodep
  for (size_type i = 0; i < num_edges; ++i)
  {
    vertex_list[f[i]]->edge_list.unsafe_push_back(edge_list[id + i]);
    if (!DIRECTION::is_directed())
    {
      vertex_list[t[i]]->edge_list.unsafe_push_back(edge_list[id + i]);
    }
  }

  mt_write(nEdges, id + num_edges);
}

/***/

template <typename DIRECTION>
void
adjacency_list<DIRECTION>::removeEdges(size_type num_edges, size_type* e)
{
  size_type id = mt_readfe(nEdges);

  // Remove the edges from the source vertices' adjacency lists.
  #pragma mta assert parallel
  for (size_type i = 0; i < num_edges; ++i)
  {
    edge_list[e[i]]->from->removeEdge(edge_list[e[i]]);
  }

  // If the graph is undirected, remove the edges from the target vertices'
  // adjacency lists.
  if (!DIRECTION::is_directed())
  {
    #pragma mta assert parallel
    for (size_type i = 0; i < num_edges; ++i)
    {
      edge_list[e[i]]->to->removeEdge(edge_list[e[i]]);
    }
  }

  // Remove the edges from the edge list that are in the area that is to be
  // deleted and delete their memory.
  #pragma mta assert nodep
  for (size_type i = 0; i < num_edges; ++i)
  {
    if (e[i] >= id - num_edges)
    {
      delete edge_list[e[i]];
      edge_list[e[i]] = 0;
    }
  }

  // Remove the edges from the edge list that are not in the area that is
  // to be deleted and delete their memory.
  size_type moved_edges = 0;
  #pragma mta assert nodep
  for (size_type i = 0; i < num_edges; ++i)
  {
    if (e[i] < id - num_edges)
    {
      // Find the next last edge in the area to be deleted that hasn't already
      // been deleted.
      size_type nme = mt_incr(moved_edges, 1);
      while (edge_list[id - 1 - nme] == 0) nme = mt_incr(moved_edges, 1);

      delete edge_list[e[i]];

      // Move the last edge into the deleted edge's position.
      edge_list[e[i]] = edge_list[id - 1 - nme];

      // Change the moved edge to have the id of the edge it replaced.
      edge_list[e[i]]->id = e[i];
    }
  }

  edge_list.resize(id - num_edges);

  mt_write(nEdges, id - num_edges);
}

/***/

template <typename DIRECTION>
template <typename T>
void
adjacency_list<DIRECTION>::removeEdges(
    size_type num_edges, Edge** e,
    edge_property_map<adjacency_list<DIRECTION>, T>* emap)
{
  typedef adjacency_list<DIRECTION> Graph;
  size_type id = mt_readfe(nEdges);

  // Remove the edges from the source vertices' adjacency lists.
  #pragma mta assert parallel
  for (size_type i = 0; i < num_edges; ++i)
  {
    e[i]->from->removeEdge(e[i]);
  }

  // If the graph is undirected, remove the edges from the target vertices'
  // adjacency lists.
  if (!DIRECTION::is_directed())
  {
    #pragma mta assert parallel
    for (size_type i = 0; i < num_edges; ++i) { 
         e[i]->to->removeEdge(e[i]); 
    }
  }

  // Remove the edges from the edge list that are in the area that is to be
  // deleted.
  #pragma mta assert nodep
  for (size_type i = 0; i < num_edges; ++i)
  {
    if (e[i]->id >= id - num_edges) edge_list[e[i]->id] = 0;
  }

  // Remove the edges from the edge list that are not in the area that is
  // to be deleted.
  size_type moved_edges = 0;
  #pragma mta assert nodep
  for (size_type i = 0; i < num_edges; ++i)
  {
    if (e[i]->id < id - num_edges)
    {
      // Find the next last edge in the area to be deleted that hasn't already
      // been deleted.
      size_type nme = mt_incr(moved_edges, 1);
      while (edge_list[id - 1 - nme] == 0) {
        nme = mt_incr(moved_edges, 1);
      }

      // Move the last edge into the deleted edge's position.
      T mval = (*emap)[detail::al_edge_adapter<Graph>(edge_list[id - 1 - nme])];
      edge_list[e[i]->id] = edge_list[id - 1 - nme];

      // Change the moved edge to have the id of the edge it replaced.
      edge_list[e[i]->id]->id = e[i]->id;

      // Bring along property_map_value.
      (*emap)[detail::al_edge_adapter<Graph>(e[i])] = mval;
    }
  }

  edge_list.resize(id - num_edges);

  mt_write(nEdges, id - num_edges);

  // Delete the memory associated with the edges.
  #pragma mta assert nodep
  for (size_type i = 0; i < num_edges; ++i) delete e[i];

}

/***/

template <typename DIRECTION>
void
adjacency_list<DIRECTION>::deep_copy(const adjacency_list& rhs)
{
  // Sets the number of vertices and edges and the edge direction.
  nVertices = rhs.vertex_list.size();
  nEdges = rhs.edge_list.size();

  vertex_list.resize(nVertices);
  edge_list.resize(nEdges);

  // Add the vertices.
  #pragma mta assert parallel
  for (size_type i = 0; i < nVertices; ++i)
  {
    vertex_list[i] = new Vertex(i);

    // Set the initial size of the vertices edge lists to the number of edges
    // they will contain.
    vertex_list[i]->edge_list.resize(rhs.vertex_list[i]->edge_list.size());
  }

  // Add the edges.
  #pragma mta assert parallel
  for (size_type i = 0; i < nEdges; ++i)
  {
    edge_list[i] = new Edge(vertex_list[rhs.edge_list[i]->from->id],
                            vertex_list[rhs.edge_list[i]->to->id], i);
  }

  // Put the edges in the vertex edge lists.
  for (size_type i = 0; i < nVertices; ++i)
  {
    for (size_type j = 0; j < rhs.vertex_list[i]->edge_list.size(); ++j)
    {
      vertex_list[i]->edge_list[j] =
        edge_list[rhs.vertex_list[i]->edge_list[j]->id];
    }
  }
}

/***/

template <typename DIRECTION>
inline
typename adjacency_list<DIRECTION>::size_type
num_vertices(const adjacency_list<DIRECTION>& g)
{
  return g.get_order();
}

/***/

template <typename DIRECTION>
inline
typename adjacency_list<DIRECTION>::size_type
num_edges(const adjacency_list<DIRECTION>& g)
{
  return g.get_size();
}

/***/

template <typename DIRECTION>
inline
typename adjacency_list<DIRECTION>::vertex_descriptor
source(const typename adjacency_list<DIRECTION>::edge_descriptor& e,
       const adjacency_list<DIRECTION>& g)
{
  return e.first;
}

/***/

template <typename DIRECTION>
inline
typename adjacency_list<DIRECTION>::vertex_descriptor
target(const typename adjacency_list<DIRECTION>::edge_descriptor& e,
       const adjacency_list<DIRECTION>& g)
{
  return e.second;
}

/***/

template <typename DIRECTION>
inline
typename adjacency_list<DIRECTION>::size_type
degree(const typename adjacency_list<DIRECTION>::vertex_descriptor& v,
       const adjacency_list<DIRECTION>& g)
{
  return v->edge_list.size();
}

/***/

template <typename DIRECTION>
inline
typename adjacency_list<DIRECTION>::size_type
out_degree(const typename
           adjacency_list<DIRECTION>::vertex_descriptor& v,
           const adjacency_list<DIRECTION>& g)
{
  return v->edge_list.size();
}

/***/

template <typename DIRECTION>
inline
typename adjacency_list<DIRECTION>::vertex_iterator
vertices(const adjacency_list<DIRECTION>& g)
{
  return g.vertices();
}

/***/

template <typename DIRECTION>
inline
typename adjacency_list<DIRECTION>::adjacency_iterator
adjacent_vertices(
    const typename adjacency_list<DIRECTION>::vertex_descriptor& v,
    const adjacency_list<DIRECTION>& g)
{
  return g.adjacent_vertices(v);
}

/***/

template <typename DIRECTION>
inline
typename adjacency_list<DIRECTION>::edge_iterator
edges(const adjacency_list<DIRECTION>& g)
{
  return g.edges();
}

/***/

template <typename DIRECTION>
inline
typename adjacency_list<DIRECTION>::out_edge_iterator
out_edges(const typename adjacency_list<DIRECTION>::vertex_descriptor& v,
          const adjacency_list<DIRECTION>& g)
{
  return g.out_edges(v);
}

/***/

template <typename DIRECTION>
inline
typename adjacency_list<DIRECTION>::thread_vertex_iterator
thread_vertices(typename adjacency_list<DIRECTION>::size_type pos,
                const adjacency_list<DIRECTION>& g)
{
  return g.thread_vertices(pos);
}

/***/

template <typename DIRECTION>
inline
typename adjacency_list<DIRECTION>::thread_adjacency_iterator
thread_adjacent_vertices(
  const typename adjacency_list<DIRECTION>::vertex_descriptor& v,
  typename adjacency_list<DIRECTION>::size_type pos,
  const adjacency_list<DIRECTION>& g)
{
  return g.thread_adjacent_vertices(v, pos);
}

/***/

template <typename DIRECTION>
inline
typename adjacency_list<DIRECTION>::thread_edge_iterator
thread_edges(typename adjacency_list<DIRECTION>::size_type pos,
             const adjacency_list<DIRECTION>& g)
{
  return g.thread_edges(pos);
}

/***/

template <typename DIRECTION>
inline
typename adjacency_list<DIRECTION>::thread_out_edge_iterator
thread_out_edges(const typename adjacency_list<DIRECTION>::vertex_descriptor& v,
                 typename adjacency_list<DIRECTION>::size_type pos,
                 const adjacency_list<DIRECTION>& g)
{
  return g.thread_out_edges(v, pos);
}

/***/

template <typename DIRECTION>
inline
typename adjacency_list<DIRECTION>::vertex_descriptor
null_vertex(const adjacency_list<DIRECTION>& g)
{
  return NULL;
}

/***/

template <typename DIRECTION>
inline
typename adjacency_list<DIRECTION>::edge_descriptor
null_edge(const adjacency_list<DIRECTION>& g)
{
  return typename adjacency_list<DIRECTION>::edge_descriptor();
}

/***/

template <typename ITERATOR, typename DIRECTION>
inline
bool
is_valid(ITERATOR& iter, typename adjacency_list<DIRECTION>::size_type p,
         const adjacency_list<DIRECTION>& tg)
{
  return true;
}

/***/

template <typename DIRECTION>
inline
bool
is_directed(const adjacency_list<DIRECTION>& g)
{
  return DIRECTION::is_directed();
}

/***/

template <typename DIRECTION>
inline
bool
is_undirected(const adjacency_list<DIRECTION>& g)
{
  return !is_directed(g);
}

/***/

template <typename DIRECTION>
inline
bool
is_bidirectional(const adjacency_list<DIRECTION>& g)
{
  return DIRECTION::is_bidirectional();
}

/***/

template <typename DIRECTION>
class vertex_id_map<adjacency_list<DIRECTION> > :
  public put_get_helper<typename adjacency_list<DIRECTION>::size_type,
                        vertex_id_map<adjacency_list<DIRECTION> > > {
public:
  typedef typename adjacency_list<DIRECTION>::vertex_descriptor key_type;
  typedef typename adjacency_list<DIRECTION>::size_type value_type;

  vertex_id_map() {}

  value_type operator[] (const key_type& k) const { return k->id; }
};

template <typename DIRECTION>
class vertex_id_map<const adjacency_list<DIRECTION> > :
  public put_get_helper<typename adjacency_list<DIRECTION>::size_type,
                        vertex_id_map<const adjacency_list<DIRECTION> > > {
public:
  typedef typename adjacency_list<DIRECTION>::vertex_descriptor key_type;
  typedef typename adjacency_list<DIRECTION>::size_type value_type;

  vertex_id_map() {}

  value_type operator[] (const key_type& k) const { return k->id; }
};


/***/

template <typename DIRECTION>
class edge_id_map<adjacency_list<DIRECTION> > :
  public put_get_helper<typename adjacency_list<DIRECTION>::size_type,
                        edge_id_map<adjacency_list<DIRECTION> > > {
public:
  typedef typename adjacency_list<DIRECTION>::edge_descriptor key_type;
  typedef typename adjacency_list<DIRECTION>::size_type value_type;

  edge_id_map() {}

  value_type operator[] (const key_type& k) const { return k.id; }
};

template <typename DIRECTION>
class edge_id_map<const adjacency_list<DIRECTION> > :
  public put_get_helper<typename adjacency_list<DIRECTION>::size_type,
                        edge_id_map<const adjacency_list<DIRECTION> > > {
public:
  typedef typename adjacency_list<DIRECTION>::edge_descriptor key_type;
  typedef typename adjacency_list<DIRECTION>::size_type value_type;

  edge_id_map() {}

  value_type operator[] (const key_type& k) const { return k.id; }
};

/***/

template <typename Graph>
struct default_hash_func<detail::al_edge_adapter<Graph> > {
  hash_size_type operator()(const detail::al_edge_adapter<Graph>& key) const
  { return integer_hash_func(key.id); }
};

/***/

template <typename DIRECTION>
inline
void init(typename adjacency_list<DIRECTION>::size_type n,
          typename adjacency_list<DIRECTION>::size_type m,
          typename adjacency_list<DIRECTION>::size_type* srcs,
          typename adjacency_list<DIRECTION>::size_type* dests,
          adjacency_list<DIRECTION>& g)
{
  return g.init(n, m, srcs, dests);
}

/***/

template <typename DIRECTION>
inline
void clear(adjacency_list<DIRECTION>& g)
{
  return g.clear();
}

/***/

template <typename DIRECTION>
inline
typename adjacency_list<DIRECTION>::vertex_descriptor
add_vertex(adjacency_list<DIRECTION>& g)
{
  return g.addVertex();
}

/***/

template <typename DIRECTION>
inline
void
add_vertices(typename adjacency_list<DIRECTION>::size_type num_verts,
             adjacency_list<DIRECTION>& g)
{
  g.addVertices(num_verts);
}

/***/

template <typename DIRECTION, typename T>
inline
void
remove_vertices(typename adjacency_list<DIRECTION>::size_type num_verts,
                typename adjacency_list<DIRECTION>::size_type* v,
                adjacency_list<DIRECTION> & g)
{
  g.removeVertices(num_verts, v);
}

/***/

template <typename DIRECTION, typename T>
inline
void
remove_vertices(typename adjacency_list<DIRECTION>::size_type num_verts,
                typename adjacency_list<DIRECTION>::size_type* v,
                edge_property_map<adjacency_list<DIRECTION>, T>& emap,
                adjacency_list<DIRECTION> & g)
{
  g.removeVertices(num_verts, v, &emap);
}

/***/

template <typename DIRECTION>
inline
void
remove_vertices(typename adjacency_list<DIRECTION>::size_type num_verts,
                typename adjacency_list<DIRECTION>::vertex_descriptor *verts,
                adjacency_list<DIRECTION> & g)
{
  g.removeVertices(num_verts, verts);
}

/***/

template <typename DIRECTION>
inline
pair<typename adjacency_list<DIRECTION>::edge_descriptor, bool>
add_edge(typename adjacency_list<DIRECTION>::vertex_descriptor f,
         typename adjacency_list<DIRECTION>::vertex_descriptor t,
         adjacency_list<DIRECTION>& g)
{
  return pair<typename adjacency_list<DIRECTION>::edge_descriptor,
              bool>(g.addEdge(f, t), true);
}

/***/

template <typename DIRECTION>
inline
pair<typename adjacency_list<DIRECTION>::edge_descriptor, bool>
add_edge(typename adjacency_list<DIRECTION>::size_type f,
         typename adjacency_list<DIRECTION>::size_type t,
         adjacency_list<DIRECTION>& g)
{
  return pair<typename adjacency_list<DIRECTION>::edge_descriptor,
              bool>(g.addEdge(f, t), true);
}

/***/

template <typename DIRECTION>
inline
void
add_edges(typename adjacency_list<DIRECTION>::size_type num_edges,
          typename adjacency_list<DIRECTION>::vertex_descriptor* f,
          typename adjacency_list<DIRECTION>::vertex_descriptor* t,
          adjacency_list<DIRECTION> & g)
{
  g.addEdges(num_edges, f, t);
}

/***/

template <typename DIRECTION>
inline
void
add_edges(typename adjacency_list<DIRECTION>::size_type num_edges,
          typename adjacency_list<DIRECTION>::size_type* f,
          typename adjacency_list<DIRECTION>::size_type* t,
          adjacency_list<DIRECTION> & g)
{
  g.addEdges(num_edges, f, t);
}

/***/

template <typename DIRECTION>
inline
void
remove_edges(typename adjacency_list<DIRECTION>::size_type num_edges,
             typename adjacency_list<DIRECTION>::size_type* e,
             adjacency_list<DIRECTION> & g)
{
  g.removeEdges(num_edges, e);
}

/***/

template <typename DIRECTION, typename T>
inline
void
remove_edges(typename adjacency_list<DIRECTION>::size_type num_edges,
             typename adjacency_list<DIRECTION>::size_type* e,
             edge_property_map<adjacency_list<DIRECTION>, T>& emap,
             adjacency_list<DIRECTION> & g)
{
  g.removeEdges(num_edges, e, &emap);
}

}

#endif
