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
/*! \file write_dimacs.hpp

    \brief Writes a graph to a dimacs formatted file.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 5/14/2010
*/
/****************************************************************************/

#ifndef MTGL_WRITE_DIMACS_HPP
#define MTGL_WRITE_DIMACS_HPP

#include <vector>
#include <sstream>

#include <mtgl/mtgl_adapter.hpp>

namespace mtgl {

namespace detail {

int write_dimacs_comp_pair(const pair<int, int>& a, const pair<int, int>& b)
{
  if (a.first == b.first)
  {
    return a.second < b.second;
  }
  else
  {
    return a.first < b.first;
  }
}

}

/// Writes a graph to a dimacs formatted file.
template <typename Graph>
bool write_dimacs(Graph& g, const char* filename)
{
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
  typedef typename graph_traits<Graph>::size_type size_type;

#ifdef DEBUG
  mt_timer timer;
  #pragma mta fence
  timer.start();
#endif

  vertex_id_map<Graph> vipm = get(_vertex_id_map, g);

  size_type order = num_vertices(g);
  size_type size = num_edges(g);

  size_type line_size = 65;
  size_type buf_size = 81 + size * line_size;
  char* buf = (char*) malloc(buf_size * sizeof(char));

  // Write the header line.
  std::ostringstream line;
  line << "p sp " << order << " " << size << std::endl;

  size_type header_size = line.str().size();
  for (size_type i = 0; i < header_size; ++i) buf[i] = line.str()[i];

  char* newbuf = buf + header_size;
  buf_size = header_size + size * line_size;

#ifdef DEBUG
  #pragma mta fence
  timer.stop();
  std::cout << "Header line write time: " << std::setw(10)
            << std::fixed << std::setprecision(6) << timer.getElapsedSeconds()
            << std::endl;

  #pragma mta fence
  timer.start();
#endif

  edge_iterator edgs = edges(g);

  // Write the edge lines.
  #pragma mta assert parallel
  for (size_type i = 0; i < size; ++i)
  {
    edge_descriptor e = edgs[i];
    vertex_descriptor src = source(e, g);
    vertex_descriptor trg = target(e, g);

    long sid = get(vipm, src);
    long tid = get(vipm, trg);

    sprintf(newbuf + line_size * i, "a %20ld %20ld                    0\n",
            sid + 1, tid + 1);
  }

#ifdef DEBUG
  #pragma mta fence
  timer.stop();
  std::cout << " Edge lines write time: " << std::setw(10)
            << std::fixed << std::setprecision(6) << timer.getElapsedSeconds()
            << std::endl;

  #pragma mta fence
  timer.start();
#endif

  bool success = write_array(filename, buf, buf_size);

#ifdef DEBUG
  #pragma mta fence
  timer.stop();
  std::cout << "       File write time: " << std::setw(10)
            << std::fixed << std::setprecision(6) << timer.getElapsedSeconds()
            << std::endl;
#endif

  free(buf);

  return success;
}

/// Writes a graph to a dimacs formatted file.
template <typename Graph, typename T>
bool write_dimacs(Graph& g, const char* filename, dynamic_array<T>& values)
{
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
  typedef typename graph_traits<Graph>::size_type size_type;

  if (values.size() == 0)
  {
    return write_dimacs(g, filename);
  }

#ifdef DEBUG
  mt_timer timer;
  #pragma mta fence
  timer.start();
#endif

  vertex_id_map<Graph> vipm = get(_vertex_id_map, g);

  size_type order = num_vertices(g);
  size_type size = num_edges(g);

  if (values.size() != size)
  {
    std::cerr << "Error in write_dimacs(): \"The values array must "
              << "contain num_edges entries.\""
              << std::endl;
    return false;
  }

  if (!std::numeric_limits<T>::is_specialized || 
      !std::numeric_limits<T>::is_integer)
  {
    std::cerr << "Error in write_dimacs(): \"The values array type "
              << "must be an integer type.\"" << std::endl;
    return false;
  }

  size_type line_size = 65;
  size_type buf_size = 81 + size * line_size;
  char* buf = (char*) malloc(buf_size * sizeof(char));

  // Write the header line.
  std::ostringstream line;
  line << "p sp " << order << " " << size << std::endl;

  size_type header_size = line.str().size();
  for (size_type i = 0; i < header_size; ++i) buf[i] = line.str()[i];

  char* newbuf = buf + header_size;
  buf_size = header_size + size * line_size;

#ifdef DEBUG
  #pragma mta fence
  timer.stop();
  std::cout << "Header line write time: " << std::setw(10)
            << std::fixed << std::setprecision(6) << timer.getElapsedSeconds()
            << std::endl;

  #pragma mta fence
  timer.start();
#endif

  edge_iterator edgs = edges(g);

  // Write the edge lines.
  #pragma mta assert parallel
  for (size_type i = 0; i < size; ++i)
  {
    edge_descriptor e = edgs[i];
    vertex_descriptor src = source(e, g);
    vertex_descriptor trg = target(e, g);

    long sid = get(vipm, src);
    long tid = get(vipm, trg);
    long val = values[i];

    sprintf(newbuf + line_size * i, "a %20ld %20ld %20ld\n",
            sid + 1, tid + 1, val);
  }

#ifdef DEBUG
  #pragma mta fence
  timer.stop();
  std::cout << " Edge lines write time: " << std::setw(10)
            << std::fixed << std::setprecision(6) << timer.getElapsedSeconds()
            << std::endl;

  #pragma mta fence
  timer.start();
#endif

  bool success = write_array(filename, buf, buf_size);

#ifdef DEBUG
  #pragma mta fence
  timer.stop();
  std::cout << "       File write time: " << std::setw(10)
            << std::fixed << std::setprecision(6) << timer.getElapsedSeconds()
            << std::endl;
#endif

  free(buf);

  return success;
}

}

#endif
