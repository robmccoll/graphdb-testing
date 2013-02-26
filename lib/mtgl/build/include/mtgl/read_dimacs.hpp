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
/*! \file read_dimacs.hpp

    \brief Parses a DIMACS graph file and creates an MTGL graph.

    \author Joe Crobak
    \author Eric Goodman (elgoodm@sandia.gov)
    \author Greg Mackey (gemacke@sandia.gov)

    \date 12/4/2007

    This function reads the dimacs graph format (.gr files) as described at
    http://www.dis.uniroma1.it/~challenge9/format.shtml.
*/
/****************************************************************************/

#ifndef MTGL_READ_DIMACS_HPP
#define MTGL_READ_DIMACS_HPP

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <iterator>

#include <mtgl/util.hpp>
#include <mtgl/algorithm.hpp>
#include <mtgl/dynamic_array.hpp>
#include <mtgl/mtgl_io.hpp>

namespace mtgl {

/// Parses a DIMACS graph file and creates an MTGL graph.
template <typename Graph, typename WGT>
bool read_dimacs(Graph& g, const char* filename, dynamic_array<WGT>& weights)
{
  typedef typename graph_traits<Graph>::size_type size_type;

#ifdef DEBUG
  mt_timer timer;
  mt_timer timer2;
  timer.start();
#endif

  long buflen;
  char* buf = read_array<char>(filename, buflen);

  if (buf == NULL) return false;

#ifdef DEBUG
  timer.stop();
  std::cout << "                 File read time: " << std::setw(10)
            << std::fixed << std::setprecision(6) << timer.getElapsedSeconds()
            << std::endl;

  timer.start();
  timer2.start();
#endif

  // Find the problem line and get the problem size.
  // Problem line expected to be in the format of:
  //    'p sp num_nodes num_edges'

  // Find the starting position in buf for the problem line.
  long num_intro_lines = 1;
  long pls = 0;

  if (buf[0] != 'p')
  {
    for ( ; pls < buflen - 1 && !(buf[pls] == '\n' && buf[pls+1] == 'p');
         ++pls)
    {
      if (buf[pls] == '\n') ++num_intro_lines;
    }

    ++num_intro_lines;
    ++pls;
  }

  if (pls == buflen)
  {
    std::cerr << std::endl << "Error: Did not find problem line." << std::endl;
    return false;
  }

  // Find the starting position for the first edge line.
  long els = pls;
  while (buf[els] != '\n') ++els;
  if (buf[els] != 'a')
  {
    for ( ; els < buflen - 1 && !(buf[els] == '\n' && buf[els+1] == 'a');
         ++els)
    {
      if (buf[els] == '\n') ++num_intro_lines;
    }

    ++num_intro_lines;
    ++els;
  }

  // Replace all the '\n' and '\r' characters with '\0' so that we have a
  // single array that is a concatenation of a bunch of propery formed
  // C strings.
  #pragma mta noalias *buf
  for (long i = 0; i < buflen; ++i)
  {
    if (buf[i] == '\n' || buf[i] == '\r') buf[i] = '\0';
  }

  size_type num_vertices;
  size_type num_edges;

  bool error = false;
  std::istringstream buf_iss(buf + pls);
  std::istream_iterator<std::string> buf_iter(buf_iss);
  if (*buf_iter != "p") error = true;
  ++buf_iter;
  if (*buf_iter != "sp") error = true;
  ++buf_iter;
  std::stringstream(*buf_iter) >> num_vertices;
  ++buf_iter;
  if (buf_iter == std::istream_iterator<std::string>()) error = true;
  std::stringstream(*buf_iter) >> num_edges;

  if (error)
  {
    std::cerr << std::endl << "Error: Problem line misformatted (expected p "
              << "sp NODES EDGES)." << std::endl;
    return false;
  }

#ifdef DEBUG
  timer2.stop();
  std::cout << "Problem line find and read time: " << std::setw(10)
            << std::fixed << std::setprecision(6) << timer2.getElapsedSeconds()
            << std::endl;

  timer2.start();
#endif

  // Find start of each edge line.
  size_type num_lines = 0;
  long* start_positions = (long*) malloc(sizeof(long) * num_edges);

  #pragma mta assert nodep
  for (long i = els - 1; i < buflen - 1; ++i)
  {
    if (buf[i] == '\0' && buf[i+1] == 'a')
    {
      size_type pos = mt_incr(num_lines, 1);
      start_positions[pos] = i + 1;
    }
  }

  if (num_lines != num_edges)
  {
    std::cerr << std::endl << "Error: Number of edges in file doesn't match "
              << "number from problem description" << std::endl
              << "       line." << std::endl;
    return false;
  }

#ifdef DEBUG
  timer2.stop();
  std::cout << "      Start positions find time: " << std::setw(10)
            << std::fixed << std::setprecision(6) << timer2.getElapsedSeconds()
            << std::endl;
#endif

#ifdef __MTA__
#ifdef DEBUG
  timer2.start();
#endif

  sort(start_positions, num_edges);

#ifdef DEBUG
  timer2.stop();
  std::cout << "      Start positions sort time: " << std::setw(10)
            << std::fixed << std::setprecision(6) << timer2.getElapsedSeconds()
            << std::endl;
#endif
#endif

#ifdef DEBUG
  timer2.start();
#endif

  // Allocate storage for the sources, destinations, and weights.
  size_type* edge_heads = (size_type*) malloc(sizeof(size_type) * num_edges);
  size_type* edge_tails = (size_type*) malloc(sizeof(size_type) * num_edges);
  weights.resize(num_edges);

  // Find all the edges.  Expected format:
  //   a SRC DST CAP
  // This will be the bulk of the lines in the graph and so is parallelized
  // where possible.  The computation of line starts is so that we can have
  // no dependencies between lines in this loop.  This, of course, implies
  // that we are being slightly inefficient for the serial case, but the
  // relative speed of the single previous pass through the buffer to find
  // the line starts appears very small for relatively large sized graphs
  // (~128MB files).

  #pragma mta assert parallel
  #pragma mta assert nodep
  for (size_type i = 0; i < num_edges; ++i)
  {
    long from;
    long to;
    long weight;

    char* a = &buf[start_positions[i] + 1];
    char* b = NULL;

    from = strtol(a, &b, 10);
    to = strtol(b, &a, 10);
    weight = strtol(a, &b, 10);

    if (a == b)
    {
      std::cerr << std::endl << "Error on line " << i + num_intro_lines + 1
                << ": Too few parameters when describing edge." << std::endl;
      exit(1);
    }

    // Dimacs vertex ids are 1-based.  We need them to be 0-based.
    --from;
    --to;

    if (from < 0 || static_cast<size_type>(from) >= num_vertices)
    {
      std::cerr << std::endl << "Error on line " << i + num_intro_lines + 1
                << ": First vertex id is invalid." << std::endl;
      exit(1);
    }
    else if (to < 0 || static_cast<size_type>(to) >= num_vertices)
    {
      std::cerr << std::endl << "Error on line " << i + num_intro_lines + 1
                << ": Second vertex id is invalid." << std::endl;
      exit(1);
    }
    else if (weight < 0)
    {
      std::cerr << std::endl << "Error on line " << i + num_intro_lines + 1
                << ": Negative weight supplied." << std::endl;
      exit(1);
    }

    edge_heads[i] = static_cast<size_type>(from);
    edge_tails[i] = static_cast<size_type>(to);
    weights[i] = static_cast<WGT>(weight);
  }

#ifdef DEBUG
  timer2.stop();
  std::cout << "                 Edge read time: " << std::setw(10)
            << std::fixed << std::setprecision(6) << timer2.getElapsedSeconds()
            << std::endl;

  timer.stop();
  std::cout << "                File parse time: " << std::setw(10)
            << std::fixed << std::setprecision(6) << timer.getElapsedSeconds()
            << std::endl;

  timer.start();
#endif

  init(num_vertices, num_edges, edge_heads, edge_tails, g);

#ifdef DEBUG
  timer.stop();
  std::cout << "            Graph creation time: " << std::setw(10)
            << std::fixed << std::setprecision(6) << timer.getElapsedSeconds()
            << std::endl;
#endif

  free(start_positions);
  free(edge_heads);
  free(edge_tails);
  free(buf);

  return true;
}

template <typename Graph>
bool read_dimacs(Graph& g, const char* filename)
{
  dynamic_array<int> values;
  return read_dimacs(g, filename, values);
}

}

#endif
