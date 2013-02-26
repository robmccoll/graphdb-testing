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
/*! \file stinger_graph.hpp

    \brief An implementation of the STINGER graph structure.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 5/18/2011
*/
/****************************************************************************/

#ifndef MTGL_STINGER_GRAPH_HPP
#define MTGL_STINGER_GRAPH_HPP

#include <iostream>
#include <iomanip>

#include <mtgl/util.hpp>

#define EDGES_IN_BLOCK 256

namespace mtgl {

template <typename stinger_graph>
class stinger_edge {
public:
  typedef typename stinger_graph::l_vertex_t l_vertex_t;
  typedef typename stinger_graph::eweight_t eweight_t;
  typedef typename stinger_graph::timestamp_t timestamp_t;

  stinger_edge() : to_id(-1) {}

  stinger_edge(l_vertex_t vid, eweight_t w, timestamp_t t1, timestamp_t t2) :
    to_id(vid), weight(w), time1(t1), time2(t2) {}

  stinger_edge(const stinger_edge& other) :
    to_id(other.to_id), weight(other.weight), time1(other.time1),
    time2(other.time2) {}

  l_vertex_t to_id;
  eweight_t weight;
  timestamp_t time1;
  timestamp_t time2;
};

template <typename stinger_graph>
class stinger_edge_block {
public:
  typedef typename stinger_graph::size_type size_type;
  typedef typename stinger_graph::l_vertex_t l_vertex_t;
  typedef typename stinger_graph::edge_t edge_t;
  typedef typename stinger_graph::etype_t etype_t;
  typedef typename stinger_graph::timestamp_t timestamp_t;

  stinger_edge_block(etype_t et, l_vertex_t vid) :
    next(0), type(et), from_id(vid), num_edges(0), smallest_t(0),
    largest_t(0) {}

  ~stinger_edge_block() {}

  stinger_edge_block* next;
  etype_t type;
  l_vertex_t from_id;
  size_type num_edges;
  timestamp_t smallest_t;
  timestamp_t largest_t;
  edge_t edges[EDGES_IN_BLOCK];
};

template <typename stinger_graph>
class stinger_vertex {
public:
  typedef typename stinger_graph::size_type size_type;
  typedef typename stinger_graph::p_vertex_t p_vertex_t;
  typedef typename stinger_graph::vtype_t vtype_t;
  typedef typename stinger_graph::vweight_t vweight_t;

  stinger_vertex() : pv_id(-1) {}

  stinger_vertex(p_vertex_t p, vtype_t t, vweight_t w) :
    pv_id(p), type(t), weight(w), in_deg(0), out_deg(0), edge_blocks(0) {}

  stinger_vertex(const stinger_vertex& other) :
    pv_id(other.pv_id), type(other.type), weight(other.weight),
    in_deg(other.in_deg), out_deg(other.out_deg),
    edge_blocks(other.edge_blocks) {}

  stinger_vertex & operator=(const stinger_vertex& other)
  {
    pv_id = other.pv_id;
    type = other.type;
    weight = other.weight;
    in_deg = other.in_deg;
    out_deg = other.out_deg;
    edge_blocks = other.edge_blocks;

    return *this;
  }

  p_vertex_t pv_id;
  vtype_t type;
  vweight_t weight;
  size_type in_deg;
  size_type out_deg;
  stinger_edge_block<stinger_graph>* edge_blocks;
};

template <typename stinger_graph>
class stinger_etype_entry {
public:
  typedef typename stinger_graph::size_type size_type;

  stinger_etype_entry() : num_blocks(0), edge_blocks(0) {}

  size_type num_blocks;
  stinger_edge_block<stinger_graph>** edge_blocks;
};


//template <typename type_t, typename weight_t, typename ts_t, typename mapper_t>
template <typename type_t, typename weight_t, typename ts_t>
class stinger_graph {
public:
  typedef int64_t size_type;
  typedef int64_t l_vertex_t;
  typedef int64_t p_vertex_t;
  typedef type_t vtype_t;
  typedef type_t etype_t;
  typedef weight_t vweight_t;
  typedef weight_t eweight_t;
  typedef ts_t timestamp_t;

  typedef stinger_vertex<stinger_graph> vertex_t;
  typedef stinger_edge<stinger_graph> edge_t;
  typedef stinger_edge_block<stinger_graph> edge_block_t;
  typedef stinger_etype_entry<stinger_graph> etype_entry_t;

//  stinger_graph(mapper_t& m, size_type net = 1) :
//    n_vertices(0), n_edges(0), n_etypes(net), vertices(0), edge_types(0),
//    initialized_by_mmap(false), mapper(m)

  stinger_graph(size_type net) : n_vertices(0), n_edges(0), n_etypes(net),
                                 vertices(0), edge_types(0), edge_blocks(0),
                                 initialized_by_mmap(false) {}

  ~stinger_graph() { clear(); }

  #pragma mta no inline
  void clear()
  {
    if (initialized_by_mmap)
    {
      initialized_by_mmap = false;
    }
    else
    {
      size_type net = n_etypes;

      // Delete the memory for the edge blocks.
      if (edge_blocks)
      {
        size_type ne = n_edges;
        #pragma mta assert parallel
        for (size_type i = 0; i < ne; ++i) edge_blocks[i].~edge_block_t();

        free(edge_blocks);
      }

      // Delete the memory for the edge type arrays.
      if (edge_types)
      {
        #pragma mta assert parallel
        for (size_type i = 0; i < net; ++i)
        {
          if (edge_types[i].edge_blocks) free(edge_types[i].edge_blocks);
        }

        #pragma mta assert par_newdelete
        delete [] edge_types;
      }

      // Delete the memory for the vertices.
      if (vertices)
      {
        size_type nv = n_vertices;
        #pragma mta assert parallel
        for (size_type i = 0; i < nv; ++i) vertices[i].~vertex_t();

        free(vertices);
      }
    }

    n_vertices = 0;
    n_edges = 0;
    vertices = 0;
    edge_types = 0;
    edge_blocks = 0;
  }

  #pragma mta no inline
  void init_vertices(l_vertex_t n, vtype_t* types, vweight_t* weights)
  {
    #pragma mta noalias *types
    #pragma mta noalias *weights

#ifdef DEBUG
    std::cout << std::endl << "Vertices: " << n << std::endl;
    for (size_type i = 0; i < n; ++i)
    {
      std::cout << "(" << std::setw(4) << i << ", " << std::setw(4) << types[i]
                << ", " << std::setw(4) << weights[i] << ")" << std::endl;
    }
#endif

    vertices = (vertex_t*) malloc(n * sizeof(vertex_t));

    // Insert the vertices.
    #pragma mta assert nodep
    for (l_vertex_t i = 0; i < n; ++i)
    {
      new(vertices + i) vertex_t(i, types[i], weights[i]);
    }

/*
    // Add the mapper entries.
    #pragma mta assert nodep
    for (l_vertex_t i = 0; i < n; ++i) mapper.insert(i, i);
*/

    n_vertices = n;
  }

  #pragma mta no inline
  void init_edges(size_type m, l_vertex_t* srcs, l_vertex_t* dests,
                  etype_t* types, eweight_t* weights)
  {
    #pragma mta noalias *srcs
    #pragma mta noalias *dests
    #pragma mta noalias *types
    #pragma mta noalias *weights

#ifdef DEBUG
    std::cout << std::endl << "Edges: " << m << std::endl;
    for (size_type i = 0; i < m; ++i)
    {
      std::cout << "(" << std::setw(4) << srcs[i] << ", " << std::setw(4)
                << dests[i] << ", " << std::setw(4) << types[i] << ", "
                << std::setw(4) << weights[i] << ")" << std::endl;
    }
#endif

    n_edges = m;
    size_type n = n_vertices;
    size_type net = n_etypes;

    // Determine the number of edge blocks of each edge type for each vertex.
    // First, have to know the number of edges of each type for each vertex.
    size_type total_vtypes = n * net;
    size_type* vertex_eblock_counts =
      (size_type*) malloc(sizeof(size_type) * total_vtypes);
    for (size_type i = 0; i < total_vtypes; ++i) vertex_eblock_counts[i] = 0;

    for (size_type i = 0; i < m; ++i)
    {
      l_vertex_t entry = srcs[i] * net + types[i];
      mt_incr(vertex_eblock_counts[entry], 1);
    }

#ifdef DEBUG
    std::cout << std::endl << "Number of edges of each type for each vertex:"
              << std::endl;
    for (size_type i = 0; i < n; ++i)
    {
      std::cout << i << ": ";
      for (size_type j = 0; j < net; ++j)
      {
        std::cout << " " << vertex_eblock_counts[i * net + j];
      }
      std::cout << std::endl;
    }
#endif

    for (size_type i = 0; i < total_vtypes; ++i)
    {
      vertex_eblock_counts[i] = (vertex_eblock_counts[i] + EDGES_IN_BLOCK - 1) /
                                EDGES_IN_BLOCK;
    }

#ifdef DEBUG
    std::cout << std::endl << "Number of edge blocks of each type for "
              << "each vertex:" << std::endl;
    for (size_type i = 0; i < n; ++i)
    {
      std::cout << i << ": ";
      for (size_type j = 0; j < net; ++j)
      {
        std::cout << " " << vertex_eblock_counts[i * net + j];
      }
      std::cout << std::endl;
    }
#endif

    // Determine the number of edge blocks of each type.
    size_type* etype_block_counts =
      (size_type*) malloc(sizeof(size_type) * net);
    for (size_type i = 0; i < net; ++i) etype_block_counts[i] = 0;

    for (size_type i = 0; i < total_vtypes; ++i)
    {
      mt_incr(etype_block_counts[i % net], vertex_eblock_counts[i]);
    }

#ifdef DEBUG
    std::cout << std::endl << "Number of edge blocks of each type:"
              << std::endl;
    for (size_type i = 0; i < net; ++i)
    {
      std::cout << i << ": " << etype_block_counts[i] << std::endl;
    }
#endif

    // Declare memory for the edge type arrays.
    #pragma mta assert par_newdelete
    edge_types = new etype_entry_t[net];

    #pragma mta assert nodep
    for (size_type i = 0; i < net; ++i)
    {
      edge_types[i].edge_blocks =
        (edge_block_t**) malloc(sizeof(edge_block_t*) * etype_block_counts[i]);
    }

    // Determine the total number of edge blocks.
    size_type total_edge_blocks = 0;
    for (size_type i = 0; i < n_etypes; ++i)
    {
      total_edge_blocks += etype_block_counts[i];
    }

    // Declare the memory for all the edge blocks.
    edge_blocks =
      (edge_block_t*) malloc(total_edge_blocks * sizeof(edge_block_t));

    total_edge_blocks = 0;

    // Create the edge blocks for each vertex.  This will parallelize the
    // outer loop over each vertex.
    #pragma mta assert nodep
    #pragma mta assert parallel
    for (size_type i = 0; i < n; ++i)
    {
      edge_block_t** cur = &(vertices[i].edge_blocks);

      for (size_type j = 0; j < net; ++j)
      {
        for (size_type k = 0; k < vertex_eblock_counts[i * net + j]; ++k)
        {
          // Add the edge block to the vertex's adjacency list.
          size_type my_pos = mt_incr(total_edge_blocks, 1);
          *cur = new(edge_blocks + my_pos) edge_block_t(j, i);

          // Add the edge block to the appropriate edge types array.
          size_type my_entry = mt_incr(edge_types[j].num_blocks, 1);
          edge_types[j].edge_blocks[my_entry] = *cur;

          cur = &((*cur)->next);
        }
      }
    }

    // Put the edges in the edge blocks.
    #pragma mta assert nodep
    for (size_type i = 0; i < m; ++i)
    {
      // Find the first open entry in an edge block of the right type in the
      // source vertex's list, to add the edge.

      // Find the first edge block of the right type.
      edge_block_t* cur = vertices[srcs[i]].edge_blocks;
      while (cur && cur->type != types[i]) cur = cur->next;
//      assert(cur);

      // Go through the edge blocks until you find an open slot.
      size_type my_slot = mt_incr(cur->num_edges, 1);
      while (my_slot >= EDGES_IN_BLOCK && cur->next)
      {
        cur = cur->next;
        my_slot = mt_incr(cur->num_edges, 1);
      }
//      assert(my_slot < EDGES_IN_BLOCK);

      // Insert the edge into the slot.
      cur->edges[my_slot].to_id = dests[i];
      cur->edges[my_slot].weight = weights[i];

      // Count this edge for the in and out degrees.
      mt_incr(vertices[srcs[i]].out_deg, 1);
      mt_incr(vertices[dests[i]].in_deg, 1);
    }

    // The num_edges field of the full edge blocks is going to be too high.
    // Go through all the edge blocks and set the num_edges field of the full
    // ones to be EDGES_IN_BLOCK.
    #pragma mta assert nodep
    for (size_type i = 0; i < net; ++i)
    {
      size_type num_blocks = edge_types[i].num_blocks;
      #pragma mta assert nodep
      for (size_type j = 0; j < num_blocks; ++j)
      {
        if (edge_types[i].edge_blocks[j]->num_edges > EDGES_IN_BLOCK)
        {
          edge_types[i].edge_blocks[j]->num_edges = EDGES_IN_BLOCK;
        }
      }
    }

    free(vertex_eblock_counts);
    free(etype_block_counts);
  }

  vertex_t& get_vertex(l_vertex_t v) const { return vertices[v]; }

  size_type num_edges() const { return n_edges; }
  size_type num_vertices() const { return n_vertices; }

  size_type out_degree(l_vertex_t v) const { return vertices[v].out_deg; }
  size_type in_degree(l_vertex_t v) const { return vertices[v].in_deg; }

//  mapper_t& get_mapper(void) const { return mapper; }

  unsigned long get_mmap_size()
  {
    // Total size:
    //
    // HEADER                 TYPE              NUM VARS
    //   mmap type            unsigned long        1
    //   mmap size            unsigned long        1
    //   order                unsigned long        1
    //   size                 unsigned long        1
    //
    // BODY                   TYPE              DIRECTED
    //   vertices             vertex_t          n_vertices
    //   edge_types           etype_entry_t     n_etypes
    //   edge_types_arrays    edge_block_t*     total_edge_blocks
    //   edge_blocks          edge_block_t      total_edge_blocks

    size_type total_edge_blocks = 0;
    for (size_type i = 0; i < n_etypes; ++i)
    {
      total_edge_blocks += edge_types[i].num_blocks;
    }

    return 4 * sizeof(unsigned long) + sizeof(uintptr_t) + sizeof(size_type) +
           n_vertices * sizeof(vertex_t) + n_etypes * sizeof(etype_entry_t) +
           total_edge_blocks * (sizeof(edge_block_t*) + sizeof(edge_block_t));
  }

  void write_mmap(void* mapped_mem, const unsigned long graph_type)
  {
    unsigned long* ul_mapped_mem = reinterpret_cast<unsigned long*>(mapped_mem);

    // Write the mmap type, mmap size, graph size, and graph order to the
    // mapped memory.
    ul_mapped_mem[0] = graph_type;
    ul_mapped_mem[1] = get_mmap_size();
    ul_mapped_mem[2] = n_vertices;
    ul_mapped_mem[3] = n_edges;

    // Write the base pointer to the memory.
    uintptr_t* base_ptr = reinterpret_cast<uintptr_t*>(ul_mapped_mem + 4);
    *base_ptr = reinterpret_cast<uintptr_t>(ul_mapped_mem);

    // Write the number of edge types to the mapped memory.
    size_type* n_etypes_ptr = reinterpret_cast<size_type*>(base_ptr + 1);
    *n_etypes_ptr = n_etypes;

    // Copy the vertices array to the mapped memory.
    vertex_t* vertices_ptr = reinterpret_cast<vertex_t*>(n_etypes_ptr + 1);
    memcpy(vertices_ptr, vertices, n_vertices * sizeof(vertex_t));

    // Set up the edge_types arrays in the mapped memory.
    etype_entry_t* edge_types_ptr =
      reinterpret_cast<etype_entry_t*>(vertices_ptr + n_vertices);

    edge_block_t** edge_types_blocks_ptr =
      reinterpret_cast<edge_block_t**>(edge_types_ptr + n_etypes);

    // Get the total number of edge blocks.
    size_type total_edge_blocks = 0;
    for (size_type i = 0; i < n_etypes; ++i)
    {
      total_edge_blocks += edge_types[i].num_blocks;
    }

    // Copy the edge blocks to the mmapped memory.
    edge_block_t* edge_blocks_ptr = reinterpret_cast<edge_block_t*>(
        edge_types_blocks_ptr + total_edge_blocks);
    memcpy(edge_blocks_ptr, edge_blocks,
           total_edge_blocks * sizeof(edge_block_t));

    for (size_type i = 0; i < n_etypes; ++i)
    {
      edge_types_ptr[i].num_blocks = edge_types[i].num_blocks;
      edge_types_ptr[i].edge_blocks = edge_types_blocks_ptr;
      edge_types_blocks_ptr += edge_types[i].num_blocks;
    }

    // Get the offset for all the edge block pointers.
    ptrdiff_t ptr_offset = reinterpret_cast<char*>(edge_blocks_ptr) -
                           reinterpret_cast<char*>(edge_blocks);

    // Add the edge block pointer offset to all the vertex adjacency list heads.
    #pragma mta assert parallel
    for (size_type i = 0; i < n_vertices; ++i)
    {
      // Only add the offset if it's a valid pointer.  Null values stay null.
      if (vertices_ptr[i].edge_blocks)
      {
        vertices_ptr[i].edge_blocks =
          reinterpret_cast<edge_block_t*>(
            reinterpret_cast<char*>(vertices[i].edge_blocks) + ptr_offset);
      }
      else
      {
        vertices_ptr[i].edge_blocks = 0;
      }
    }

    // Add the edge block pointer offset to all the pointers in the edge
    // block type arrays.
    #pragma mta assert parallel
    #pragma mta assert nodep
    for (size_type i = 0; i < n_etypes; ++i)
    {
      #pragma mta assert parallel
      #pragma mta assert nodep
      for (size_type j = 0; j < edge_types_ptr[i].num_blocks; ++j)
      {
        edge_types_ptr[i].edge_blocks[j] =
          reinterpret_cast<edge_block_t*>(
            reinterpret_cast<char*>(edge_types[i].edge_blocks[j]) + ptr_offset);
      }
    }

    // Add the edge block pointer offset to all the vertex adjacency list
    // heads.
    #pragma mta assert parallel
    for (size_type i = 0; i < total_edge_blocks; ++i)
    {
      // Only add the offset if it's a valid pointer.  Null values stay null.
      if (edge_blocks[i].next)
      {
        edge_blocks_ptr[i].next = 
          reinterpret_cast<edge_block_t*>(
            reinterpret_cast<char*>(edge_blocks[i].next) + ptr_offset);
      }
      else
      {
        edge_blocks_ptr[i].next = 0;
      }
    }
  }

  void read_mmap(void* mapped_mem)
  {
    clear();

    unsigned long* ul_mapped_mem = reinterpret_cast<unsigned long*>(mapped_mem);

    // Set the size and order of the graph.
    n_vertices = ul_mapped_mem[2];
    n_edges = ul_mapped_mem[3];

    // Get the offset for all the pointers.
    uintptr_t* base_ptr = reinterpret_cast<uintptr_t*>(ul_mapped_mem + 4);
    ptrdiff_t ptr_offset = reinterpret_cast<char*>(ul_mapped_mem) -
                           reinterpret_cast<char*>(*base_ptr);

    // Set the number of edge types for the graph.
    size_type* n_etypes_ptr = reinterpret_cast<size_type*>(base_ptr + 1);
    n_etypes = *n_etypes_ptr;

    // Get pointers to the locations in the mapped memory for the arrays.
    vertex_t* vertices_ptr = reinterpret_cast<vertex_t*>(n_etypes_ptr + 1);
    etype_entry_t* edge_types_ptr =
      reinterpret_cast<etype_entry_t*>(vertices_ptr + n_vertices);

    // Update the pointers in the edge types arrays with the offset.
    // NOTE: This modifies data in the mmapped region of memory.
    #pragma mta assert parallel
    for (size_type i = 0; i < n_etypes; ++i)
    {
      edge_types_ptr[i].edge_blocks =
        reinterpret_cast<edge_block_t**>(
          reinterpret_cast<char*>(edge_types_ptr[i].edge_blocks) + ptr_offset);
    }

    #pragma mta assert nodep
    #pragma mta assert parallel
    for (size_type i = 0; i < n_etypes; ++i)
    {
      #pragma mta assert nodep
      #pragma mta assert parallel
      for (size_type j = 0; j < edge_types_ptr[i].num_blocks; ++j)
      {
        edge_types_ptr[i].edge_blocks[j] =
          reinterpret_cast<edge_block_t*>(
            reinterpret_cast<char*>(edge_types_ptr[i].edge_blocks[j]) +
                                    ptr_offset);
      }
    }

    // Update the edge block pointers in the vertices adjacency lists.  This
    // updates the pointers in the edge blocks, themselves, too.
    // NOTE: This modifies data in the mmapped region of memory.
    #pragma mta assert nodep
    #pragma mta assert parallel
    for (size_type i = 0; i < n_vertices; ++i)
    {
      edge_block_t** cur = &(vertices_ptr[i].edge_blocks);

      while (*cur)
      {
        *cur = reinterpret_cast<edge_block_t*>(
                 reinterpret_cast<char*>(*cur) + ptr_offset);

        cur = &((*cur)->next);
      }
    }

    // Set the pointers to the arrays for the graph to point to the mapped
    // memory.
    vertices = vertices_ptr;
    edge_types = edge_types_ptr;

    // We don't really care about the edge_blocks pointer in a graph
    // initialized by mmapped memory.  It's only used to allocate and
    // deallocate memory in a single chunk.
    edge_blocks = 0;

    // Write my pointer to the mapped mem as the new base pointer;
    *base_ptr = reinterpret_cast<uintptr_t>(ul_mapped_mem);

    initialized_by_mmap = true;
  }

private:
  size_type n_vertices;
  size_type n_edges;
  size_type n_etypes;

  vertex_t* vertices;
  etype_entry_t* edge_types;
  edge_block_t* edge_blocks;
  bool initialized_by_mmap;
//  mapper_t& mapper;
};

}

#undef EDGES_IN_BLOCK

#endif
