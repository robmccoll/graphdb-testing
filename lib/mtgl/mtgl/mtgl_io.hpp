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
/*! \file mtgl_io.hpp

    \brief This include file defines or includes all the necessary functions
           for I/O using the MTGL I/O interface.

    \author Greg Mackey (gemacke@sandia.gov)
    \author Eric Goodman (elgoodm@sandia.gov)

    \date 2/11/2010
*/
/****************************************************************************/

#ifndef MTGL_MTGL_IO_HPP
#define MTGL_MTGL_IO_HPP

#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <cstdio>

#include <mtgl/snap_util.h>

namespace mtgl {

/*! \brief A convenience method that given a file, it returns the size of
           the file.  Returns 0 if something goes wrong.  Assumes that
           snap_init has already been called.  snap_util.h should hopefully
           make this compatible with non XMT systems.

    \param file The file that will be checked.

    \author Eric Goodman (elgoodm@sandia.gov)
    \date March 2010
*/
inline
uint64_t get_file_size(const char *file)
{
  luc_error_t err;
  int64_t snap_error = 0;

  // Try to us any fsworker
  luc_endpoint_id_t endpoint = SNAP_ANY_SW;


  // Used to store the results of the stat call.
  snap_stat_buf stat_buf;
  err = snap_stat(const_cast<char*>(file), endpoint, &stat_buf, &snap_error);

  if (SNAP_ERR_OK != err)
  {
    //Try to do things the old way with SW_WORKER_EP
    char* fsw = getenv("SWORKER_EP");
    if (fsw != NULL) endpoint = strtoul(fsw, NULL, 16);
    err = snap_stat(const_cast<char*>(file), endpoint, &stat_buf, &snap_error);

#ifdef DEBUG
    std::cerr << "endpoint_string: " << fsw << std::endl
              << "       endpoint: " << std::hex << std::setw(16)
              << std::setfill('0') << endpoint << std::endl;
#endif

    if (SNAP_ERR_OK != err)
    {
      std::cerr << "Failed to retrieve info on file " << file << ". Error "
                << err << "." << std::endl;
      return 0;
    }
  }

  return stat_buf.st_size;
}

/*! \brief Reads a binary file and puts the contents of the file into a
           buffer of type array_t.

    \param filename The name of the file to be read in via snap_restore
    \param size The number of array_t elements in the returned array.  It
                is set to 0 if there are problems.

    \returns Returns an array of array_t elements.

    \author Eric Goodman (elgoodm@sandia.gov)
    \date March 2010
*/
template <typename array_t, typename size_type>
inline
array_t* read_array(const char* filename, size_type& array_size)
{
  array_t* array;
  array_size = 0;

  // Initializing snapshot libraries.
  if (snap_init() != SNAP_ERR_OK)
  {
    perror("Can't initialize libsnapshot.\n");
    return NULL;
  }

  // Getting the size of the file.
  uint64_t file_size = get_file_size(filename);

  if (file_size == 0) return NULL;

  // Creating the buffer to store the data from the file.
  void* buffer = (void*) malloc(file_size * sizeof(char));

  if (buffer == NULL)
  {
    perror("Failed to malloc snapshot buffer.\n");
    return NULL;
  }

  if (snap_restore(const_cast<char*>(filename),
                   buffer, file_size, NULL) != SNAP_ERR_OK)
  {
    perror("Couldn't snap_restore file.\n");
    free(buffer);
    return NULL;
  }

  // Setting the size of the array.
  array_size = file_size / sizeof(array_t);

  // Casting the buffer to be the appropriate type.
  array = (array_t*) buffer;

  return array;
}

/*! \brief Writes a binary file with the contents of the given array.

    \param filename The name of the file to be written to via snap_snapshot
    \param array The array to written
    \param size The number of array_t elements in array.

    \author Eric Goodman (elgoodm@sandia.gov)

    \date March 2010
*/
template <typename array_t, typename size_type>
inline
bool write_array(const char* filename, array_t* array, size_type array_size)
{
  // Initializing snapshot libraries.
  if (snap_init() != SNAP_ERR_OK)
  {
    std::cerr << "Error in write_array(): \"Can't initialize libsnapshot.\""
              << std::endl;
    return false;
  }

  uint64_t file_size = array_size * sizeof(array_t);
  if (snap_snapshot(const_cast<char*>(filename),
                    array, file_size, NULL) != SNAP_ERR_OK)
  {
    std::cerr << "Error in write_array(): \"Couldn't snap_snapshot file.\""
              << std::endl;
    return false;
  }

  return true;
}

/*! \brief This function performs byte swapping to switch between little and
           big endian and vice-versa on a single memory location.

    \param buf Pointer to the memory to byte swap.
    \param len Number of bytes in buf.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 2/17/2011

    As an example, let's perform byte_swapping on val, an int.

    int val = ...;

    byte_swap(&val, sizeof(int));
*/
void byte_swap(void* buf, size_t len)
{
  char* c = static_cast<char*>(buf);

  #pragma mta loop serial
  for (size_t i = 0; i < len / 2; ++i)
  {
    char t = c[i];
    c[i] = c[len - i - 1];
    c[len - i - 1] = t;
  }
}

/*! \brief This function performs byte swapping to switch between little and
           big endian and vice-versa on all the entries in an array.

    \param buf Pointer to the array to byte swap.
    \param array_len Length of buf.
    \param len Number of bytes in each entry of buf.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 9/14/2011

    As an example, let's perform byte_swapping on vals, an array of ints.

    int* vals = ...;
    int vals_length = number of integers in the array;

    byte_swap(vals, vals_length, sizeof(int));
*/
template <typename T1, typename T2>
void byte_swap(T1* buf, T2 array_len, size_t len)
{
  #pragma mta noalias *buf

  for (T2 i = 0; i < array_len; ++i)
  {
    char* c = reinterpret_cast<char*>(buf + i);

    #pragma mta loop serial
    for (size_t i = 0; i < len / 2; ++i)
    {
      char t = c[i];
      c[i] = c[len - i - 1];
      c[len - i - 1] = t;
    }
  }
}

}

#include <mtgl/read_binary.hpp>
#include <mtgl/write_binary.hpp>
#include <mtgl/read_dimacs.hpp>
#include <mtgl/write_dimacs.hpp>
#include <mtgl/read_matrix_market.hpp>
#include <mtgl/write_matrix_market.hpp>
#include <mtgl/read_mmap.hpp>
#include <mtgl/write_mmap.hpp>

#endif
