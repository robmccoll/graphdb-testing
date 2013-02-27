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
/*! \file read_binary.hpp

    \brief Reads in a graph from the srcs / dests format.

    \author Greg Mackey (gemacke@sandia.gov)
    \author Jon Berry (jberry@sandia.gov)

    \date 5/14/2010

    This code was pulled from the old restore() member function of graph
    adapters.
*/
/****************************************************************************/

#ifndef MTGL_READ_BINARY_HPP
#define MTGL_READ_BINARY_HPP

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>

#include <mtgl/snap_util.h>
#include <mtgl/mtgl_adapter.hpp>

namespace mtgl {

/// \brief Restore edge source and destination data from snapshot files.
///        The type of the records in the file is assumed to be size_type.
template <typename Graph>
bool read_binary(Graph& g, const char* src_filename, const char* dest_filename)
{
  typedef typename graph_traits<Graph>::size_type size_type;

#ifdef DEBUG
  mt_timer timer;
  timer.start();
#endif

  if (snap_init() != SNAP_ERR_OK)
  {
    perror("Can't initialize libsnapshot.\n");
    return false;
  }

  // fstat both files.
  char* swEP_string = getenv("SWORKER_EP");
  luc_endpoint_id_t swEP = 0;
  if (swEP_string != NULL) swEP = strtoul(swEP_string, NULL, 16);

#ifdef DEBUG
  if (swEP_string)
  {
    std::cout << "swEP_string: " << swEP_string << std::endl;
  }
  else
  {
    std::cout << "swEP_string: " << "NULL" << std::endl;
  }
  std::cout << "       swEP: " << std::setfill('0') << std::hex
            << std::setw(16) << swEP << std::endl << std::setfill(' ')
            << std::dec;
#endif

  snap_stat_buf statBuf;
  int64_t sn_err = 0;

  // Get size of data in src and dest files; compute number of records in
  // each file.
  if (!src_filename ||
      snap_stat(const_cast<char*>(src_filename),
                swEP, &statBuf, &sn_err) != SNAP_ERR_OK)
  {
    perror("Can't stat sources file.\n");
    return false;
  }

  size_t src_size = statBuf.st_size / sizeof(size_type);

  if (!dest_filename ||
      snap_stat (const_cast<char*>(dest_filename),
                 swEP, &statBuf, &sn_err) != SNAP_ERR_OK)
  {
    perror("Can't stat destinations file.\n");
    return false;
  }

  size_t dest_size = statBuf.st_size / sizeof(size_type);

  if (src_size != dest_size)
  {
    perror("The number of entries in the srcs and dests files are "
           "not equal.\n");
    return false;
  }

  // Allocate arrays for srcs and dests.
  size_type* srcs = (size_type*) malloc(src_size * sizeof(size_type));
  size_type* dests = (size_type*) malloc(dest_size * sizeof(size_type));

  // Restore srcs into buffer.
  if (snap_restore(const_cast<char*>(src_filename),
                   srcs, statBuf.st_size, NULL) != SNAP_ERR_OK)
  {
    perror("Couldn't snap_restore sources file.\n");
    return false;
  }

  // Restore dests into buffer.
  if (snap_restore(const_cast<char*>(dest_filename),
                   dests, statBuf.st_size, NULL) != SNAP_ERR_OK)
  {
    perror("couldn't snap_restore destinations file.\n");
    return false;
  }

#ifdef DEBUG
  timer.stop();
  std::cout << "      File read time: " << std::setw(10)
            << std::fixed << std::setprecision(6) << timer.getElapsedSeconds()
            << std::endl;
  timer.start();
#endif

  // Scan for max vertex to find the number of verices in the graph.
  size_type max_vertex = 0;
  #pragma mta block schedule
  for (size_t i = 0; i < src_size; ++i)
  {
    if (srcs[i] > max_vertex) max_vertex = srcs[i];
  }

  #pragma mta block schedule
  for (size_t i = 0; i < dest_size; ++i)
  {
    if (dests[i] > max_vertex) max_vertex = dests[i];
  }

#ifdef DEBUG
  timer.stop();
  std::cout << "Max vertex find time: " << std::setw(10)
            << std::fixed << std::setprecision(6) << timer.getElapsedSeconds()
            << std::endl;
#endif

  size_type n = max_vertex + 1;
  size_type m = src_size;

#ifdef DEBUG
  timer.start();
#endif

  init(n, m, srcs, dests, g);

#ifdef DEBUG
  timer.stop();
  std::cout << " Graph creation time: " << std::setw(10)
            << std::fixed << std::setprecision(6) << timer.getElapsedSeconds()
            << std::endl;
#endif

  free(srcs);
  free(dests);

  return true;
}

/// \brief Restore edge source and destination data from snapshot files.
///        The type of the records in the files is int_type.
template <typename Graph, typename int_type>
bool read_binary(Graph& g, const char* src_filename, const char* dest_filename,
                 int_type val)
{
  typedef typename graph_traits<Graph>::size_type size_type;

#ifdef DEBUG
  mt_timer timer;
  timer.start();
#endif

  if (snap_init() != SNAP_ERR_OK)
  {
    perror("Can't initialize libsnapshot.\n");
    return false;
  }

  // fstat both files.
  char* swEP_string = getenv("SWORKER_EP");
  luc_endpoint_id_t swEP = 0;
  if (swEP_string != NULL) swEP = strtoul(swEP_string, NULL, 16);

#ifdef DEBUG
  if (swEP_string)
  {
    std::cout << "swEP_string: " << swEP_string << std::endl;
  }
  else
  {
    std::cout << "swEP_string: " << "NULL" << std::endl;
  }
  std::cout << "       swEP: " << std::setfill('0') << std::hex
            << std::setw(16) << swEP << std::endl << std::setfill(' ')
            << std::dec;
#endif

  snap_stat_buf statBuf;
  int64_t sn_err = 0;

  // Get size of data in src and dest files; compute number of records in
  // each file.
  if (!src_filename ||
      snap_stat(const_cast<char*>(src_filename),
                swEP, &statBuf, &sn_err) != SNAP_ERR_OK)
  {
    perror("Can't stat sources file.\n");
    return false;
  }

  size_t src_size = statBuf.st_size / sizeof(int_type);

  if (!dest_filename ||
      snap_stat (const_cast<char*>(dest_filename),
                 swEP, &statBuf, &sn_err) != SNAP_ERR_OK)
  {
    perror("Can't stat destinations file.\n");
    return false;
  }

  size_t dest_size = statBuf.st_size / sizeof(int_type);

  if (src_size != dest_size)
  {
    perror("The number of entries in the srcs and dests files are "
           "not equal.\n");
    return false;
  }

  // Allocate arrays for srcs and dests.
  size_type* srcs = (size_type*) malloc(src_size * sizeof(size_type));
  size_type* dests = (size_type*) malloc(dest_size * sizeof(size_type));

  // Allocate a tmp buffer for use by restore.  (Need this tmp buffer
  // since sizeof(size_type) differs from sizeof(int_type).)
  int_type* tmpbuf = (int_type*) malloc(statBuf.st_size);

  if (sizeof(int_type) > sizeof(size_type))
  {
    perror("snap_snapshot:  File record size is greater than "
           "sizeof(size_type); the data read may be incorrect.");
  }

  // Restore srcs into buffer.
  if (snap_restore(const_cast<char*>(src_filename),
                   tmpbuf, statBuf.st_size, NULL) != SNAP_ERR_OK)
  {
    perror("Couldn't snap_restore sources file.\n");
    return false;
  }

  // Copy restored srcs into srcs array, casting to correct type.
  for (size_type i = 0; i < src_size; ++i)
  {
    srcs[i] = static_cast<size_type>(tmpbuf[i]);
  }

  // Restore dests into buffer.
  if (snap_restore(const_cast<char*>(dest_filename),
                   tmpbuf, statBuf.st_size, NULL) != SNAP_ERR_OK)
  {
    perror("couldn't snap_restore destinations file.\n");
    return false;
  }

  // Copy restored dests into dests array, casting to correct type.
  for (size_type i = 0; i < dest_size; ++i)
  {
    dests[i] = static_cast<size_type>(tmpbuf[i]);
  }

  free(tmpbuf);

#ifdef DEBUG
  timer.stop();
  std::cout << "      File read time: " << std::setw(10)
            << std::fixed << std::setprecision(6) << timer.getElapsedSeconds()
            << std::endl;
  timer.start();
#endif

  // Scan for max vertex to find the number of verices in the graph.
  size_type max_vertex = 0;
  for (size_type i = 0; i < src_size; ++i)
  {
    if (srcs[i] > max_vertex) max_vertex = srcs[i];
  }

  for (size_type i = 0; i < dest_size; ++i)
  {
    if (dests[i] > max_vertex) max_vertex = dests[i];
  }

#ifdef DEBUG
  timer.stop();
  std::cout << "Max vertex find time: " << std::setw(10)
            << std::fixed << std::setprecision(6) << timer.getElapsedSeconds()
            << std::endl;
#endif

  size_type n = max_vertex + 1;
  size_type m = src_size;

#ifdef DEBUG
  timer.start();
#endif

  init(n, m, srcs, dests, g);

#ifdef DEBUG
  timer.stop();
  std::cout << " Graph creation time: " << std::setw(10)
            << std::fixed << std::setprecision(6) << timer.getElapsedSeconds()
            << std::endl;
#endif

  free(srcs);
  free(dests);

  return true;
}

}

#endif
