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
/*! \file read_matrix_market.hpp

    \brief Parses a matrix market graph file and creates an MTGL graph.

    \author Karen Devine (kddevin@sandia.gov)
    \author Greg Mackey (gemacke@sandia.gov)

    \date 2/26/2008
*/
/****************************************************************************/

#ifndef MTGL_READ_MATRIX_MARKET_HPP
#define MTGL_READ_MATRIX_MARKET_HPP

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <iterator>

#include <mtgl/util.hpp>
#include <mtgl/mtgl_io.hpp>
#include <mtgl/algorithm.hpp>
#include <mtgl/dynamic_array.hpp>

// The XMT implementation of strtod() had horrible performance in
// read_matrix_market() (where horrible means read_matrix_market() and not
// just the bad version of strtod() ran 10x slower), so I added a version
// that didn't suck.
#ifdef __MTA__
  #include <mtgl/strtod.hpp>
#endif

namespace mtgl {

#ifdef USING_QTHREADS
namespace detail {

class replace_line_enders {
public:
  replace_line_enders(char* b) : buf(b) {}

  void operator()(const size_t start, const size_t stop)
  {
    for (size_t i = start; i != stop; ++i)
    {
      if (buf[i] == '\n' || buf[i] == '\r') buf[i] = '\0';
    }
  }

private:
  char* buf;
};

template <typename size_type>
class find_starts {
public:
  find_starts(char* b, size_type& nl, long* sp) :
    buf(b), num_lines(nl), start_positions(sp) {}

  void operator()(const size_t start, const size_t stop)
  {
    size_type local_positions[64];
    unsigned int local_num_lines = 0;

    for (size_t i = start; i < stop; ++i)
    {
      if (buf[i] == '\0' && buf[i+1] != '\0' && buf[i+1] != '%')
      {
        local_positions[local_num_lines++] = i + 1;

        if (local_num_lines == 64)
        {
          size_type pos = mt_incr(num_lines, 64);

          for (int j = 0; j < 64; ++j)
          {
            start_positions[pos+j] = local_positions[j];
          }

          local_num_lines = 0;
        }
      }
    }

    if (local_num_lines > 0)
    {
      size_type pos = mt_incr(num_lines, local_num_lines);

      for (unsigned int j = 0; j < local_num_lines; ++j)
      {
        start_positions[pos+j] = local_positions[j];
      }
    }
  }

private:
  char* buf;
  size_type& num_lines;
  long* start_positions;
};

// weight_cat:
//   0 - pattern
//   1 - integer
//   2 - real
//   3 - complex
template <typename size_type, typename T, typename WT, int weight_cat>
class parse_edges {
public:
  parse_edges(char* b, size_type nil, long* sp, size_type* eh, size_type *et,
              dynamic_array<T>& v, size_type nv, size_type ne) :
    buf(b), num_intro_lines(nil), start_positions(sp), edge_heads(eh),
    edge_tails(et), values(v), num_vertices(nv), num_edges(ne) {}

  void operator()(const size_t start, const size_t stop)
  {
    for (size_t i = start; i != stop; ++i)
    {
      long from;
      long to;
      WT value;
      WT imag;

      char* a = &buf[start_positions[i]];
      char* b = NULL;

      from = strtol(a, &b, 10);
      to = strtol(b, &a, 10);

      if (weight_cat == 1)
      {
        value = strtol(a, &b, 10);
      }
      else if (weight_cat == 2)
      {
        value = strtod(a, &b);
      }
      else if (weight_cat == 3)
      {
        value = strtod(a, &b);
        imag = strtod(b, &a);
      }

      if (a == b)
      {
        std::cerr << std::endl << "Error on line " << i + num_intro_lines + 1
                  << ": Too few parameters when describing edge."
                  << std::endl;
        exit(1);
      }

      // Matrix Market vertex ids are 1-based.  We need them to be 0-based.
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

      edge_heads[i] = static_cast<size_type>(from);
      edge_tails[i] = static_cast<size_type>(to);

      if (weight_cat != 0) values[i] = static_cast<T>(value);
      if (weight_cat == 3) values[i + num_edges] = static_cast<T>(imag);
    }
  }

private:
  char* buf;
  size_type num_intro_lines;
  long* start_positions;
  size_type* edge_heads;
  size_type* edge_tails;
  dynamic_array<T>& values;
  size_type num_vertices;
  size_type num_edges;
};

}
#endif

/*! \brief Parses a matrix market graph file and creates an MTGL graph.

    \author Karen Devine (kddevin@sandia.gov)
    \author Greg Mackey (gemacke@sandia.gov)
*/
template <typename Graph, typename T>
bool read_matrix_market(Graph& g, const char* filename,
                        dynamic_array<T>& values)
{
  typedef typename graph_traits<Graph>::size_type size_type;

#ifdef DEBUG
  mt_timer timer;
  mt_timer timer2;
  #pragma mta fence
  timer.start();
#endif

  long buflen;
  char* buf = read_array<char>(filename, buflen);

  if (buf == NULL) return false;

#ifdef DEBUG
  #pragma mta fence
  timer.stop();
  std::cout << "                 File read time: " << std::setw(10)
            << std::fixed << std::setprecision(6) << timer.getElapsedSeconds()
            << std::endl;

  #pragma mta fence
  timer.start();
  #pragma mta fence
  timer2.start();
#endif

  // Find the beginning of the problem line.
  long num_intro_lines = 0;
  long pls = 0;

  while (buf[pls] == '%' || buf[pls] == '\n')
  {
    while (buf[pls] != '\n') ++pls;    // Skip the line.

    ++num_intro_lines;
    ++pls;                            // Move to next line
  }

  // Replace all the '\n' and '\r' characters with '\0' so that we have a
  // single array that is a concatenation of a bunch of propery formed
  // C strings.
#ifdef USING_QTHREADS
  detail::replace_line_enders rle(buf);
  qt_loop_balance(0, buflen, rle);
#else
  for (long i = 0; i < buflen; ++i)
  {
    if (buf[i] == '\n' || buf[i] == '\r') buf[i] = '\0';
  }
#endif

  char mtx[256] = {0};
  char crd[256];
  char entry_type[256];
  char symmetry_format[256];

  int ret = 0;
  ret = sscanf(buf, "%%%%MatrixMarket %s %s %s %s", mtx, crd, entry_type,
               symmetry_format);

  if (ret != 4)
  {
    std::cerr << std::endl << "Error: Header line misformatted. (ret = "
              << ret << ")" << std::endl;
    return false;
  }

  // Convert to lower case.
  for (char* p = mtx; *p != '\0'; *p = tolower(*p), ++p);
  for (char* p = crd; *p != '\0'; *p = tolower(*p), ++p);
  for (char* p = entry_type; *p != '\0'; *p = tolower(*p), ++p);
  for (char* p = symmetry_format; *p != '\0'; *p = tolower(*p), ++p);

  if (strcmp(mtx, "matrix") != 0)
  {
    std::cerr << std::endl << "Error: Can only read 'matrix' objects."
              << std::endl;
    return false;
  }

  if (strcmp(crd, "coordinate") != 0)
  {
    std::cerr << std::endl << "Error: Can only read 'coordinate' format."
              << std::endl;
    return false;
  }

  if (strcmp(entry_type, "real") != 0 &&
      strcmp(entry_type, "complex") != 0 &&
      strcmp(entry_type, "integer") != 0 &&
      strcmp(entry_type, "pattern") != 0)
  {
    std::cerr << std::endl << "Error: Can only read 'real', 'complex', "
              << "'integer', or 'pattern' matrix entry types." << std::endl;
    return false;
  }

  if (strcmp(symmetry_format, "general") != 0 &&
      strcmp(symmetry_format, "symmetric") != 0 &&
      strcmp(symmetry_format, "skew-symmetric") != 0 &&
      strcmp(symmetry_format, "hermitian") != 0)
  {
    std::cerr << std::endl << "Error: Can only read 'general', 'symmetric', "
              << "'skew-symmetric', and 'hermitian' symmetry" << std::endl
              << "       formats." << std::endl;
    return false;
  }

  if (strcmp(entry_type, "pattern") == 0 &&
      strcmp(symmetry_format, "skew-symmetric") == 0)
  {
    std::cerr << std::endl << "Error: Can't have a pattern skew-symmetric "
              << "graph." << std::endl;
    return false;
  }

  if (strcmp(symmetry_format, "hermitian") == 0 &&
      strcmp(entry_type, "complex") != 0)
  {
    std::cerr << std::endl << "Error: Hermitian graphs must have a 'complex' "
              << "field type." << std::endl;
    return false;
  }

  // Read size of sparse matrix.

  size_type num_vertices;
  size_type ncol;
  size_type num_entries;

  std::istringstream buf_iss(buf + pls);
  std::istream_iterator<std::string> buf_iter(buf_iss);
  std::stringstream(*buf_iter) >> num_vertices;
  ++buf_iter;
  std::stringstream(*buf_iter) >> ncol;
  ++buf_iter;
  bool error = false;
  if (buf_iter == std::istream_iterator<std::string>()) error = true;
  std::stringstream(*buf_iter) >> num_entries;

  if (error)
  {
    std::cerr << std::endl << "Error: Problem line misformatted (expected "
              << "<num_rows> <num_cols> <num_nonzeros>)." << std::endl;
    return false;
  }

  if (num_vertices != ncol)
  {
    std::cerr << std::endl << "Error: The matrix must be square "
              << num_vertices << " x " << ncol << std::endl;
    return false;
  }

  // Move pls past the problem line.
  while (buf[pls] != '\0') ++pls;
  ++num_intro_lines;

#ifdef DEBUG
  #pragma mta fence
  timer2.stop();
  std::cout << "Problem line find and read time: " << std::setw(10)
            << std::fixed << std::setprecision(6) << timer2.getElapsedSeconds()
            << std::endl;

  #pragma mta fence
  timer2.start();
#endif

  // Find start of each edge line.
  size_type num_lines = 0;
  long* start_positions = (long*) malloc(sizeof(long) * num_entries);

#ifdef USING_QTHREADS
  detail::find_starts<size_type> fs(buf, num_lines, start_positions);
  qt_loop_balance(pls, buflen - 1, fs);
#else
  #pragma mta assert nodep
  for (long i = pls; i < buflen - 1; ++i)
  {
    if (buf[i] == '\0' && buf[i+1] != '\0' && buf[i+1] != '%')
    {
      size_type pos = mt_incr(num_lines, 1);
      start_positions[pos] = i + 1;
    }
  }
#endif

  if (num_lines != num_entries)
  {
    std::cerr << std::endl << "Error: Number of edges in file doesn't match "
              << "number from problem description" << std::endl
              << "       line." << std::endl;
    return false;
  }

#ifdef DEBUG
  #pragma mta fence
  timer2.stop();
  std::cout << "      Start positions find time: " << std::setw(10)
            << std::fixed << std::setprecision(6) << timer2.getElapsedSeconds()
            << std::endl;
#endif

#ifdef __MTA__
#ifdef DEBUG
  #pragma mta fence
  timer2.start();
#endif

  sort(start_positions, num_entries);

#ifdef DEBUG
  #pragma mta fence
  timer2.stop();
  std::cout << "      Start positions sort time: " << std::setw(10)
            << std::fixed << std::setprecision(6) << timer2.getElapsedSeconds()
            << std::endl;
#endif
#endif

#ifdef DEBUG
  #pragma mta fence
  timer2.start();
#endif

  // Find the number of edges.
  size_type num_edges = 0;
  size_type num_diags = 0;
  if (strcmp(symmetry_format, "general") == 0)
  {
    // For general symmetry format, there are num_entries edges.
    num_edges = num_entries;
  }
  else if (strcmp(symmetry_format, "skew-symmetric") == 0)
  {
    // For skew-symmetric symmetry format, there are 2 * num_entries edges.
    num_edges = 2 * num_entries;
  }
  else
  {
    //  For symmetric and Hermitian symmetry formats there are
    // num_diags + 2 * num_offdiags edges, so we need to count the number of
    // entries on the diagonal to calculate the number of edges.
    for (size_type i = 0; i < num_entries; ++i)
    {
      long from;
      long to;

      char* a = &buf[start_positions[i]];
      char* b = NULL;

      from = strtol(a, &b, 10);
      to = strtol(b, &a, 10);

      // No need to check for errors in reading from and to because strtol()
      // always returns a valid integer.  Errors will be caught later.
      num_diags += from == to;
    }

    num_edges = num_diags + 2 * (num_entries - num_diags);
  }

  // Allocate storage for the sources, destinations, and values.
  size_type* edge_heads = (size_type*) malloc(sizeof(size_type) * num_edges);
  size_type* edge_tails = (size_type*) malloc(sizeof(size_type) * num_edges);

  // Complex entry type has 2 * num_edges values to store the real and
  // imaginary parts of the complex number.  Real and integer entry types have
  // num_edges values.  Pattern entry type has no values.
  if (strcmp(entry_type, "complex") == 0)
  {
    values.resize(2 * num_edges);
  }
  else if (strcmp(entry_type, "pattern") != 0)
  {
    values.resize(num_edges);
  }

  // Read edges.

  if (strcmp(symmetry_format, "general") == 0)
  {
#ifdef USING_QTHREADS
    if (strcmp(entry_type, "pattern") == 0)
    {
      detail::parse_edges<size_type, T, long, 0>
        pe(buf, num_intro_lines, start_positions, edge_heads, edge_tails,
           values, num_vertices, num_edges);
      qt_loop_balance(0, num_entries, pe);
    }
    else if (strcmp(entry_type, "integer") == 0)
    {
      detail::parse_edges<size_type, T, long, 1>
        pe(buf, num_intro_lines, start_positions, edge_heads, edge_tails,
           values, num_vertices, num_edges);
      qt_loop_balance(0, num_entries, pe);
    }
    else if (strcmp(entry_type, "real") == 0)
    {
      detail::parse_edges<size_type, T, double, 2>
        pe(buf, num_intro_lines, start_positions, edge_heads, edge_tails,
           values, num_vertices, num_edges);
      qt_loop_balance(0, num_entries, pe);
    }
    else if (strcmp(entry_type, "complex") == 0)
    {
      detail::parse_edges<size_type, T, double, 3>
        pe(buf, num_intro_lines, start_positions, edge_heads, edge_tails,
           values, num_vertices, num_edges);
      qt_loop_balance(0, num_entries, pe);
    }
#else
    if (strcmp(entry_type, "real") == 0)
    {
      #pragma mta assert parallel
      for (size_type i = 0; i < num_entries; ++i)
      {
        long from;
        long to;
        double value;

        char* a = &buf[start_positions[i]];
        char* b = NULL;

        from = strtol(a, &b, 10);
        to = strtol(b, &a, 10);
        value = strtod(a, &b);

        if (a == b)
        {
          std::cerr << std::endl << "Error on line " << i + num_intro_lines + 1
                    << ": Too few parameters when describing edge."
                    << std::endl;
          exit(1);
        }

        // Matrix Market vertex ids are 1-based.  We need them to be 0-based.
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

        edge_heads[i] = static_cast<size_type>(from);
        edge_tails[i] = static_cast<size_type>(to);
        values[i] = static_cast<T>(value);
      }
    }
    else if (strcmp(entry_type, "complex") == 0)
    {
      #pragma mta assert parallel
      for (size_type i = 0; i < num_entries; ++i)
      {
        long from;
        long to;
        double real;
        double imag;

        char* a = &buf[start_positions[i]];
        char* b = NULL;

        from = strtol(a, &b, 10);
        to = strtol(b, &a, 10);
        real = strtod(a, &b);
        imag = strtod(b, &a);

        if (a == b)
        {
          std::cerr << std::endl << "Error on line " << i + num_intro_lines + 1
                    << ": Too few parameters when describing edge."
                    << std::endl;
          exit(1);
        }

        // Matrix Market vertex ids are 1-based.  We need them to be 0-based.
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

        edge_heads[i] = static_cast<size_type>(from);
        edge_tails[i] = static_cast<size_type>(to);
        values[i] = static_cast<T>(real);
        values[i + num_edges] = static_cast<T>(imag);
      }
    }
    else if (strcmp(entry_type, "integer") == 0)
    {
      #pragma mta assert parallel
      for (size_type i = 0; i < num_entries; ++i)
      {
        long from;
        long to;
        long value;

        char* a = &buf[start_positions[i]];
        char* b = NULL;

        from = strtol(a, &b, 10);
        to = strtol(b, &a, 10);
        value = strtol(a, &b, 10);

        if (a == b)
        {
          std::cerr << std::endl << "Error on line " << i + num_intro_lines + 1
                    << ": Too few parameters when describing edge."
                    << std::endl;
          exit(1);
        }

        // Matrix Market vertex ids are 1-based.  We need them to be 0-based.
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

        edge_heads[i] = static_cast<size_type>(from);
        edge_tails[i] = static_cast<size_type>(to);
        values[i] = static_cast<T>(value);
      }
    }
    else if (strcmp(entry_type, "pattern") == 0)
    {
      #pragma mta assert parallel
      for (size_type i = 0; i < num_entries; ++i)
      {
        long from;
        long to;

        char* a = &buf[start_positions[i]];
        char* b = NULL;

        from = strtol(a, &b, 10);
        to = strtol(b, &a, 10);

        if (a == b)
        {
          std::cerr << std::endl << "Error on line " << i + num_intro_lines + 1
                    << ": Too few parameters when describing edge."
                    << std::endl;
          exit(1);
        }

        // Matrix Market vertex ids are 1-based.  We need them to be 0-based.
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

        edge_heads[i] = static_cast<size_type>(from);
        edge_tails[i] = static_cast<size_type>(to);
      }
    }
#endif
  }
  else if (strcmp(symmetry_format, "symmetric") == 0)
  {
    if (strcmp(entry_type, "real") == 0)
    {
      long my_num_diags = 0;
      long my_num_offdiags = 0;

      #pragma mta assert parallel
      for (size_type i = 0; i < num_entries; ++i)
      {
        long from;
        long to;
        double value;

        char* a = &buf[start_positions[i]];
        char* b = NULL;

        from = strtol(a, &b, 10);
        to = strtol(b, &a, 10);
        value = strtod(a, &b);

        if (a == b)
        {
          std::cerr << std::endl << "Error on line " << i + num_intro_lines + 1
                    << ": Too few parameters when describing edge."
                    << std::endl;
          exit(1);
        }

        // Matrix Market vertex ids are 1-based.  We need them to be 0-based.
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

        if (to == from)
        {
          // Diagonal edges only get one entry in the matrix.
          long pos = mt_incr(my_num_diags, 1);
          edge_heads[pos] = static_cast<size_type>(from);
          edge_tails[pos] = static_cast<size_type>(to);
          values[pos] = static_cast<T>(value);
        }
        else
        {
          // Off diagonal edges get two entries in the matrix.
          long pos = mt_incr(my_num_offdiags, 1);
          edge_heads[num_diags + 2 * pos] = static_cast<size_type>(from);
          edge_tails[num_diags + 2 * pos] = static_cast<size_type>(to);
          values[num_diags + 2 * pos] = static_cast<T>(value);
          edge_heads[num_diags + 2 * pos + 1] = static_cast<size_type>(to);
          edge_tails[num_diags + 2 * pos + 1] = static_cast<size_type>(from);
          values[num_diags + 2 * pos + 1] = static_cast<T>(value);
        }
      }
    }
    else if (strcmp(entry_type, "complex") == 0)
    {
      long my_num_diags = 0;
      long my_num_offdiags = 0;

      #pragma mta assert parallel
      for (size_type i = 0; i < num_entries; ++i)
      {
        long from;
        long to;
        double real;
        double imag;

        char* a = &buf[start_positions[i]];
        char* b = NULL;

        from = strtol(a, &b, 10);
        to = strtol(b, &a, 10);
        real = strtod(a, &b);
        imag = strtod(b, &a);

        if (a == b)
        {
          std::cerr << std::endl << "Error on line " << i + num_intro_lines + 1
                    << ": Too few parameters when describing edge."
                    << std::endl;
          exit(1);
        }

        // Matrix Market vertex ids are 1-based.  We need them to be 0-based.
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

        if (to == from)
        {
          // Diagonal edges only get one entry in the matrix.
          long pos = mt_incr(my_num_diags, 1);
          edge_heads[pos] = static_cast<size_type>(from);
          edge_tails[pos] = static_cast<size_type>(to);
          values[pos] = static_cast<T>(real);
          values[pos + num_edges] = static_cast<T>(imag);
        }
        else
        {
          // Off diagonal edges get two entries in the matrix.
          long pos = mt_incr(my_num_offdiags, 1);
          edge_heads[num_diags + 2 * pos] = static_cast<size_type>(from);
          edge_tails[num_diags + 2 * pos] = static_cast<size_type>(to);
          values[num_diags + 2 * pos] = static_cast<T>(real);
          values[num_diags + 2 * pos + num_edges] = static_cast<T>(imag);
          edge_heads[num_diags + 2 * pos + 1] = static_cast<size_type>(to);
          edge_tails[num_diags + 2 * pos + 1] = static_cast<size_type>(from);
          values[num_diags + 2 * pos + 1] = static_cast<T>(real);
          values[num_diags + 2 * pos + 1 + num_edges] = static_cast<T>(imag);
        }
      }
    }
    else if (strcmp(entry_type, "integer") == 0)
    {
      long my_num_diags = 0;
      long my_num_offdiags = 0;

      #pragma mta assert parallel
      for (size_type i = 0; i < num_entries; ++i)
      {
        long from;
        long to;
        long value;

        char* a = &buf[start_positions[i]];
        char* b = NULL;

        from = strtol(a, &b, 10);
        to = strtol(b, &a, 10);
        value = strtol(a, &b, 10);

        if (a == b)
        {
          std::cerr << std::endl << "Error on line " << i + num_intro_lines + 1
                    << ": Too few parameters when describing edge."
                    << std::endl;
          exit(1);
        }

        // Matrix Market vertex ids are 1-based.  We need them to be 0-based.
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

        if (to == from)
        {
          // Diagonal edges only get one entry in the matrix.
          long pos = mt_incr(my_num_diags, 1);
          edge_heads[pos] = static_cast<size_type>(from);
          edge_tails[pos] = static_cast<size_type>(to);
          values[pos] = static_cast<T>(value);
        }
        else
        {
          // Off diagonal edges get two entries in the matrix.
          long pos = mt_incr(my_num_offdiags, 1);
          edge_heads[num_diags + 2 * pos] = static_cast<size_type>(from);
          edge_tails[num_diags + 2 * pos] = static_cast<size_type>(to);
          values[num_diags + 2 * pos] = static_cast<T>(value);
          edge_heads[num_diags + 2 * pos + 1] = static_cast<size_type>(to);
          edge_tails[num_diags + 2 * pos + 1] = static_cast<size_type>(from);
          values[num_diags + 2 * pos + 1] = static_cast<T>(value);
        }
      }
    }
    else if (strcmp(entry_type, "pattern") == 0)
    {
      long my_num_diags = 0;
      long my_num_offdiags = 0;

      #pragma mta assert parallel
      for (size_type i = 0; i < num_entries; ++i)
      {
        long from;
        long to;

        char* a = &buf[start_positions[i]];
        char* b = NULL;

        from = strtol(a, &b, 10);
        to = strtol(b, &a, 10);

        if (a == b)
        {
          std::cerr << std::endl << "Error on line " << i + num_intro_lines + 1
                    << ": Too few parameters when describing edge."
                    << std::endl;
          exit(1);
        }

        // Matrix Market vertex ids are 1-based.  We need them to be 0-based.
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

        if (to == from)
        {
          // Diagonal edges only get one entry in the matrix.
          long pos = mt_incr(my_num_diags, 1);
          edge_heads[pos] = static_cast<size_type>(from);
          edge_tails[pos] = static_cast<size_type>(to);
        }
        else
        {
          // Off diagonal edges get two entries in the matrix.
          long pos = mt_incr(my_num_offdiags, 1);
          edge_heads[num_diags + 2 * pos] = static_cast<size_type>(from);
          edge_tails[num_diags + 2 * pos] = static_cast<size_type>(to);
          edge_heads[num_diags + 2 * pos + 1] = static_cast<size_type>(to);
          edge_tails[num_diags + 2 * pos + 1] = static_cast<size_type>(from);
        }
      }
    }
  }
  else if (strcmp(symmetry_format, "skew-symmetric") == 0)
  {
    if (strcmp(entry_type, "real") == 0)
    {
      #pragma mta assert parallel
      for (size_type i = 0; i < num_entries; ++i)
      {
        long from;
        long to;
        double value;

        char* a = &buf[start_positions[i]];
        char* b = NULL;

        from = strtol(a, &b, 10);
        to = strtol(b, &a, 10);
        value = strtod(a, &b);

        if (a == b)
        {
          std::cerr << std::endl << "Error on line " << i + num_intro_lines + 1
                    << ": Too few parameters when describing edge."
                    << std::endl;
          exit(1);
        }

        // Matrix Market vertex ids are 1-based.  We need them to be 0-based.
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
        else if (from == to)
        {
          std::cerr << std::endl << "Error on line " << i + num_intro_lines + 1
                    << ": Can't have diagonal entry in skew-symmetric matrix."
                    << std::endl;
          exit(1);
        }

        edge_heads[2 * i] = static_cast<size_type>(from);
        edge_tails[2 * i] = static_cast<size_type>(to);
        values[2 * i] = static_cast<T>(value);

        edge_heads[2 * i + 1] = static_cast<size_type>(to);
        edge_tails[2 * i + 1] = static_cast<size_type>(from);
        values[2 * i + 1] = static_cast<T>(-value);
      }
    }
    else if (strcmp(entry_type, "complex") == 0)
    {
      #pragma mta assert parallel
      for (size_type i = 0; i < num_entries; ++i)
      {
        long from;
        long to;
        double real;
        double imag;

        char* a = &buf[start_positions[i]];
        char* b = NULL;

        from = strtol(a, &b, 10);
        to = strtol(b, &a, 10);
        real = strtod(a, &b);
        imag = strtod(b, &a);

        if (a == b)
        {
          std::cerr << std::endl << "Error on line " << i + num_intro_lines + 1
                    << ": Too few parameters when describing edge."
                    << std::endl;
          exit(1);
        }

        // Matrix Market vertex ids are 1-based.  We need them to be 0-based.
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
        else if (from == to)
        {
          std::cerr << std::endl << "Error on line " << i + num_intro_lines + 1
                    << ": Can't have diagonal entry in skew-symmetric matrix."
                    << std::endl;
          exit(1);
        }

        edge_heads[2 * i] = static_cast<size_type>(from);
        edge_tails[2 * i] = static_cast<size_type>(to);
        values[2 * i] = static_cast<T>(real);
        values[2 * i + num_edges] = static_cast<T>(imag);

        edge_heads[2 * i + 1] = static_cast<size_type>(to);
        edge_tails[2 * i + 1] = static_cast<size_type>(from);
        values[2 * i + 1] = static_cast<T>(-real);
        values[2 * i + 1 + num_edges] = static_cast<T>(-imag);
      }
    }
    else if (strcmp(entry_type, "integer") == 0)
    {
      #pragma mta assert parallel
      for (size_type i = 0; i < num_entries; ++i)
      {
        long from;
        long to;
        long value;

        char* a = &buf[start_positions[i]];
        char* b = NULL;

        from = strtol(a, &b, 10);
        to = strtol(b, &a, 10);
        value = strtol(a, &b, 10);

        if (a == b)
        {
          std::cerr << std::endl << "Error on line " << i + num_intro_lines + 1
                    << ": Too few parameters when describing edge."
                    << std::endl;
          exit(1);
        }

        // Matrix Market vertex ids are 1-based.  We need them to be 0-based.
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
        else if (from == to)
        {
          std::cerr << std::endl << "Error on line " << i + num_intro_lines + 1
                    << ": Can't have diagonal entry in skew-symmetric matrix."
                    << std::endl;
          exit(1);
        }

        edge_heads[2 * i] = static_cast<size_type>(from);
        edge_tails[2 * i] = static_cast<size_type>(to);
        values[2 * i] = static_cast<T>(value);

        edge_heads[2 * i + 1] = static_cast<size_type>(to);
        edge_tails[2 * i + 1] = static_cast<size_type>(from);
        values[2 * i + 1] = static_cast<T>(-value);
      }
    }
  }
  else if (strcmp(symmetry_format, "hermitian") == 0)
  {
    // No need to check entry_type because a hermitian graph must have a
    // complex entry_type.

    long my_num_diags = 0;
    long my_num_offdiags = 0;

    #pragma mta assert parallel
    for (size_type i = 0; i < num_entries; ++i)
    {
      long from;
      long to;
      double real;
      double imag;

      char* a = &buf[start_positions[i]];
      char* b = NULL;

      from = strtol(a, &b, 10);
      to = strtol(b, &a, 10);
      real = strtod(a, &b);
      imag = strtod(b, &a);

      if (a == b)
      {
        std::cerr << std::endl << "Error on line " << i + num_intro_lines + 1
                  << ": Too few parameters when describing edge." << std::endl;
        exit(1);
      }

      // Matrix Market vertex ids are 1-based.  We need them to be 0-based.
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

      if (to == from)
      {
        // Diagonal edges only get one entry in the matrix.
        long pos = mt_incr(my_num_diags, 1);
        edge_heads[pos] = static_cast<size_type>(from);
        edge_tails[pos] = static_cast<size_type>(to);
        values[pos] = static_cast<T>(real);
        values[pos + num_edges] = static_cast<T>(imag);
      }
      else
      {
        // Off diagonal edges get two entries in the matrix.
        long pos = mt_incr(my_num_offdiags, 1);
        edge_heads[num_diags + 2 * pos] = static_cast<size_type>(from);
        edge_tails[num_diags + 2 * pos] = static_cast<size_type>(to);
        values[num_diags + 2 * pos] = static_cast<T>(real);
        values[num_diags + 2 * pos + num_edges] = static_cast<T>(imag);
        edge_heads[num_diags + 2 * pos + 1] = static_cast<size_type>(to);
        edge_tails[num_diags + 2 * pos + 1] = static_cast<size_type>(from);
        values[num_diags + 2 * pos + 1] = static_cast<T>(real);
        values[num_diags + 2 * pos + 1 + num_edges] = static_cast<T>(-imag);
      }
    }
  }

#ifdef DEBUG
  #pragma mta fence
  timer2.stop();
  std::cout << "                 Edge read time: " << std::setw(10)
            << std::fixed << std::setprecision(6) << timer2.getElapsedSeconds()
            << std::endl;

  #pragma mta fence
  timer.stop();
  std::cout << "                File parse time: " << std::setw(10)
            << std::fixed << std::setprecision(6) << timer.getElapsedSeconds()
            << std::endl;

  #pragma mta fence
  timer.start();
#endif

  init(num_vertices, num_edges, edge_heads, edge_tails, g);

#ifdef DEBUG
  #pragma mta fence
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

/*! \brief Parses a matrix market graph file and creates an MTGL graph.

    \author Karen Devine (kddevin@sandia.gov)
    \author Greg Mackey (gemacke@sandia.gov)
*/
template <typename Graph>
bool read_matrix_market(Graph& g, const char* filename)
{
  dynamic_array<int> values;
  return read_matrix_market(g, filename, values);
}

}

#endif
