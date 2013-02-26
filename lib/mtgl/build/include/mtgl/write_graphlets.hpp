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
/*! \file write_graphlets.hpp

    \brief Writes a graph to a file in David Gleich's graphlets format.

    \author Jon Berry (jberry@sandia.gov)

    \date 4/14/2011
*/
/****************************************************************************/

#ifndef MTGL_WRITE_GRAPHLETS_HPP
#define MTGL_WRITE_GRAPHLETS_HPP

#include <cstdio>
#include <sstream>
#include <iomanip>
#include <limits>
#include <complex>

#include <mtgl/mtgl_adapter.hpp>

namespace mtgl {

/// Writes a graph to a file in David Gleich's graphlets format.
template <typename Graph>
bool write_graphlets(Graph& ga, const char* filename)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

  vertex_id_map<Graph> vipm = get(_vertex_id_map, ga);
  edge_id_map<Graph> eipm = get(_edge_id_map, ga);

  size_type order = num_vertices(ga);
  size_type size = num_edges(ga);

  size_type line_size = 44;   // JWB DEBUG!  was 42
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

  // Write the edge lines.
  edge_iterator edgs = edges(ga);
  #pragma mta assert parallel
  for (size_type i = 0; i < size; ++i)
  {
    edge_descriptor e = edgs[i];
    vertex_descriptor src = source(e, ga);
    vertex_descriptor trg = target(e, ga);

    size_type sid = get(vipm, src);
    size_type tid = get(vipm, trg);

    // Matrix-market files are one-based.
    std::ostringstream edge_line;

    //edge_line << std::setw(20) << sid + 1 << " " << std::setw(20) << tid + 1
    //          << std::endl;
    // JWB DEBUG: TESTING WITH GRAPHLETS
    edge_line << std::setw(20) << sid << " " << std::setw(20) << tid
              << " 1" << std::endl;

    for (int j = 0; j < line_size; ++j)
    {
      newbuf[line_size * i + j] = edge_line.str()[j];
    }
  }

  bool success = write_array(filename, buf, buf_size);

  free(buf);

  return success;
}

/// Writes a graph to a file in David Gleich's graphlets format.
template <typename Graph, typename T>
bool write_graphlets(Graph& ga, const char* filename, dynamic_array<T>& values)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

  vertex_id_map<Graph> vipm = get(_vertex_id_map, ga);
  edge_id_map<Graph> eipm = get(_edge_id_map, ga);

  size_type order = num_vertices(ga);
  size_type size = num_edges(ga);

  if (values.size() != size && values.size() != 2 * size)
  {
    std::cerr << "Error in write_graphlets(): \"The values array must "
              << "contain either num_edges or 2 * num_edges entries.\""
              << std::endl;
    return false;
  }

  if (!std::numeric_limits<T>::is_specialized)
  {
    std::cerr << "Error in write_graphlets(): \"The values array type "
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

  if (entry_type != "complex")
  {
    // This is either an integer or real type.
    // Write the edge lines.
    edge_iterator edgs = edges(ga);
    #pragma mta assert parallel
    for (size_type i = 0; i < size; ++i)
    {
      edge_descriptor e = edgs[i];
      vertex_descriptor src = source(e, ga);
      vertex_descriptor trg = target(e, ga);

      size_type sid = get(vipm, src);
      size_type tid = get(vipm, trg);

      // Matrix-market files are one-based.
      std::ostringstream edge_line;
      edge_line << std::setw(20) << sid + 1 << " " << std::setw(20) << tid + 1
                << " " << std::setw(20) << std::setprecision(12)
                << std::scientific << values[i] << std::endl;

      for (int j = 0; j < line_size; ++j)
      {
        newbuf[line_size * i + j] = edge_line.str()[j];
      }
    }
  }
  else
  {
    // This is a complex type.
    // Write the edge lines.
    edge_iterator edgs = edges(ga);
    #pragma mta assert parallel
    for (size_type i = 0; i < size; ++i)
    {
      edge_descriptor e = edgs[i];
      vertex_descriptor src = source(e, ga);
      vertex_descriptor trg = target(e, ga);

      size_type sid = get(vipm, src);
      size_type tid = get(vipm, trg);

      // Matrix-market files are one-based.
      std::ostringstream edge_line;
      edge_line << std::setw(20) << sid + 1 << " " << std::setw(20) << tid + 1
                << " " << std::setprecision(12) << std::scientific
                << std::setw(20) << values[i] << " " << std::setw(20)
                << values[size + i] << std::endl;

      for (int j = 0; j < line_size; ++j)
      {
        newbuf[line_size * i + j] = edge_line.str()[j];
      }
    }
  }

  bool success = write_array(filename, buf, buf_size);

  free(buf);

  return success;
}

}

#endif
