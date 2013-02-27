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
/*! \file write_matrix_market.hpp

    \brief Writes a graph to a matrix-market formatted file.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 5/14/2010

    We only write 'coordinate' (sparse) graphs in the 'general' symmetry
    format (aka no using formats that take advantage of symmetry).  There are
    two forms of the function.  The first form uses the 'pattern' matrix entry
    type to specify that there are no values associated with the edges.  The
    first argument gives the graph to write, and the second argument gives the
    filename to write to.

      template <typename Graph, typename T>
      write_matrix_market(Graph& g, const char* filename)

    The second form writes a value associated with the edge.  This value can
    be either integer, real, or complex.  This version has the same two first
    arguments.  The third argument is a dyanmic_array of the values to write.
    If the size of the dynamic_array is twice the number of edges, the value
    type is set to be "complex".  Otherwise, if the template type of the
    dynamic_array is an integer type, the value type is set to be "integer".
    Otherwise, the value type is set to be "real".

      template <typename Graph, typename T>
      write_matrix_market(Graph& g, const char* filename,
                          dynamic_array<T>& values)
*/
/****************************************************************************/

#ifndef MTGL_WRITE_MATRIX_MARKET_HPP
#define MTGL_WRITE_MATRIX_MARKET_HPP

#include <sstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <string>

#include <mtgl/mtgl_adapter.hpp>
#include <mtgl/util.hpp>

namespace mtgl {

/// Writes a graph to a matrix-market formatted file with no value for the
/// matrix entry.
template <typename Graph>
bool write_matrix_market(Graph& g, const char* filename)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

#ifdef DEBUG
  mt_timer timer;
  #pragma mta fence
  timer.start();
#endif

  vertex_id_map<Graph> vipm = get(_vertex_id_map, g);

  size_type order = num_vertices(g);
  size_type size = num_edges(g);

  size_type line_size = 42;
  size_type buf_size = 162 + size * line_size;
  char* buf = (char*) malloc(buf_size * sizeof(char));

  // Write the header line.
  std::ostringstream line;
  line << "%%MatrixMarket matrix coordinate pattern general" << std::endl
       << order << " " << order << " " << size << std::endl;

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
  #pragma mta assert nodep
  for (size_type i = 0; i < size; ++i)
  {
    edge_descriptor e = edgs[i];
    vertex_descriptor src = source(e, g);
    vertex_descriptor trg = target(e, g);

    long sid = get(vipm, src);
    long tid = get(vipm, trg);

    sprintf(newbuf + line_size * i, "%20ld %20ld\n", sid + 1, tid + 1);
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

/// Writes a graph to a matrix-market formatted file with a value for the
/// matrix entry.
template <typename Graph, typename T>
bool write_matrix_market(Graph& g, const char* filename,
                         dynamic_array<T>& values)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

  if (values.size() == 0)
  {
    return write_matrix_market(g, filename);
  }

#ifdef DEBUG
  mt_timer timer;
  #pragma mta fence
  timer.start();
#endif

  vertex_id_map<Graph> vipm = get(_vertex_id_map, g);

  size_type order = num_vertices(g);
  size_type size = num_edges(g);

  if (values.size() != size && values.size() != 2 * size)
  {
    std::cerr << "Error in write_matrix_market(): \"The values array must "
              << "contain either num_edges or 2 * num_edges entries.\""
              << std::endl;
    return false;
  }

  if (!std::numeric_limits<T>::is_specialized)
  {
    std::cerr << "Error in write_matrix_market(): \"The values array type "
              << "must be a built-in numerical type.\"" << std::endl;
    return false;
  }

  std::string entry_type = "real";
  size_type line_size = 63;
  if (values.size() == 2 * size)
  {
    entry_type = "complex";
    line_size = 84;
  }
  else if (std::numeric_limits<T>::is_integer)
  {
    entry_type = "integer";
  }

  size_type buf_size = 162 + size * line_size;
  char* buf = (char*) malloc(buf_size * sizeof(char));

  // Write the header line.
  std::ostringstream line;
  line << "%%MatrixMarket matrix coordinate " << entry_type.c_str()
       << " general" << std::endl
       << order << " " << order << " " << size << std::endl;

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

  if (entry_type == "complex")
  {
    // This is a complex type.

    // Write the edge lines.
    #pragma mta assert parallel
    for (size_type i = 0; i < size; ++i)
    {
      edge_descriptor e = edgs[i];
      vertex_descriptor src = source(e, g);
      vertex_descriptor trg = target(e, g);

      long sid = get(vipm, src);
      long tid = get(vipm, trg);
      double val1 = values[i];
      double val2 = values[size + i];

      sprintf(newbuf + line_size * i, "%20ld %20ld %20.12e %20.12e\n",
              sid + 1, tid + 1, val1, val2);
    }
  }
  else if (entry_type == "real")
  {
    // This is a real type.

    // Write the edge lines.
    #pragma mta assert parallel
    for (size_type i = 0; i < size; ++i)
    {
      edge_descriptor e = edgs[i];
      vertex_descriptor src = source(e, g);
      vertex_descriptor trg = target(e, g);

      long sid = get(vipm, src);
      long tid = get(vipm, trg);
      double val = values[i];

      sprintf(newbuf + line_size * i, "%20ld %20ld %20.12e\n",
              sid + 1, tid + 1, val);
    }
  }
  else
  {
    // This is an integer type.

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

      sprintf(newbuf + line_size * i, "%20ld %20ld %20ld\n",
              sid + 1, tid + 1, val);
    }
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
