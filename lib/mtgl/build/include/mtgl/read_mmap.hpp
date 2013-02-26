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
/*! \file read_mmap.hpp

    \brief Initializes a graph from an mmapped memory region.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 11/11/2010
*/
/****************************************************************************/

#ifndef MTGL_READ_MMAP_HPP
#define MTGL_READ_MMAP_HPP

#ifndef _WIN32

#include <limits>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>

#include <mtgl/mtgl_adapter.hpp>
#include <mtgl/mtgl_io.hpp>
#include <mtgl/xmt_hash_table.hpp>

namespace mtgl {

template <typename Graph>
bool read_mmap(Graph& g, const char* filename)
{
  int mmap_prots = PROT_READ | PROT_WRITE;
  int mmap_flags = MAP_SHARED;

#ifdef __MTA__
  mmap_flags |= MAP_ANON;
  int fd = open(filename, O_RDWR, S_IRUSR | S_IWUSR);
#else
  int fd = shm_open(filename, O_RDWR, S_IRUSR | S_IWUSR);
#endif

  if (fd < 0)
  {
    std::cout << "Error opening file: "<< filename << std::endl;
    return false;
  }

  unsigned long mapped_size = 4 * sizeof(unsigned long);

  // Map the memory for the header.
  void* mapped_mem = mmap(0, mapped_size, mmap_prots, mmap_flags, fd, 0);

  if (mapped_mem == MAP_FAILED)
  {
    std::cout << "Error mapping memory: "<< filename << std::endl;
    close(fd);
    return false;
  }

  unsigned long* ul_mapped_mem = reinterpret_cast<unsigned long*>(mapped_mem);

  // Read the mmap type and size from the mapped memory.
  unsigned long mmap_type = ul_mapped_mem[0];
  unsigned long mmap_size = ul_mapped_mem[1];

  // Unmap the memory, so we can map the whole thing now that we have the
  // size.
  munmap((caddr_t) mapped_mem, mapped_size);

  // Make sure the graph I'm reading has the same type as me.
  if (mmap_type != mmap_traits<Graph>::type)
  {
    std::cout << "Error, wrong mmap type: "<< filename << std::endl;
    close(fd);
    return false;
  }

  mapped_mem = mmap(0, mmap_size, mmap_prots, mmap_flags, fd, 0);

  if (mapped_mem == MAP_FAILED)
  {
    std::cout << "Error mapping memory: "<< filename << std::endl;
    close(fd);
    return false;
  }

  g.read_mmap(mapped_mem);

  close(fd);

  return true;
}

template <typename T>
unsigned long read_mmap(T*& v, const char* filename)
{
  int mmap_prots = PROT_READ;
  int mmap_flags = MAP_SHARED;

#ifdef __MTA__
  mmap_flags |= MAP_ANON;
  int fd = open(filename, O_RDONLY);
#else
  int fd = shm_open(filename, O_RDONLY, S_IRUSR);
#endif

  if (fd < 0)
  {
    std::cout << "Error opening file: "<< filename << std::endl;
    return false;
  }

  unsigned long mapped_size = 4 * sizeof(unsigned long);

  // Map the memory for the header.
  void* mapped_mem = mmap(0, mapped_size, mmap_prots, mmap_flags, fd, 0);

  if (mapped_mem == MAP_FAILED)
  {
    std::cout << "Error mapping memory: "<< filename << std::endl;
    close(fd);
    return 0;
  }

  unsigned long* ul_mapped_mem = reinterpret_cast<unsigned long*>(mapped_mem);

  // Read the mmap size and type from the mapped memory.
  unsigned long mmap_type = ul_mapped_mem[0];
  unsigned long mmap_size = ul_mapped_mem[1];
  unsigned long size = ul_mapped_mem[2];

  // Unmap the memory, so we can map the whole thing now that we have the
  // size.
  munmap((caddr_t) mapped_mem, mapped_size);

  // Make sure the vector I'm reading is of the correct type.
  if (mmap_type != mmap_traits<T>::type)
  {
    std::cout << "Error, wrong mmap type: "<< filename << std::endl;
    close(fd);
    return 0;
  }

  mapped_mem = mmap(0, mmap_size, mmap_prots, mmap_flags, fd, 0);

  if (mapped_mem == MAP_FAILED)
  {
    std::cout << "Error mapping memory: "<< filename << std::endl;
    close(fd);
    return 0;
  }

  // Set v to be the mmapped array, incrementing past the header.
  ul_mapped_mem = reinterpret_cast<unsigned long*>(mapped_mem);
  v = reinterpret_cast<T*>(ul_mapped_mem + 4);

  close(fd);

  return size;
}

template <typename K, typename T, typename HF, typename EQF>
bool read_mmap(xmt_hash_table<K, T, HF, EQF>& ht, const char* filename)
{
  int mmap_prots = PROT_READ;
  int mmap_flags = MAP_SHARED;

#ifdef __MTA__
  mmap_flags |= MAP_ANON;
  int fd = open(filename, O_RDONLY);
#else
  int fd = shm_open(filename, O_RDONLY, S_IRUSR);
#endif

  if (fd < 0)
  {
    std::cout << "Error opening file: "<< filename << std::endl;
    return false;
  }

  unsigned long mapped_size = 4 * sizeof(unsigned long);

  // Map the memory for the header.
  void* mapped_mem = mmap(0, mapped_size, mmap_prots, mmap_flags, fd, 0);

  if (mapped_mem == MAP_FAILED)
  {
    std::cout << "Error mapping memory: "<< filename << std::endl;
    close(fd);
    return false;
  }

  unsigned long* ul_mapped_mem = reinterpret_cast<unsigned long*>(mapped_mem);

  // Read the mmap type and size from the mapped memory.
  unsigned long mmap_type = ul_mapped_mem[0];
  unsigned long mmap_size = ul_mapped_mem[1];

  // Unmap the memory, so we can map the whole thing now that we have the
  // size.
  munmap((caddr_t) mapped_mem, mapped_size);

  // Make sure the graph I'm reading has the same type as me.
  if (mmap_type != mmap_traits<xmt_hash_table<K, T, HF, EQF> >::type)
  {
    std::cout << "Error, wrong mmap type: "<< filename << std::endl;
    close(fd);
    return false;
  }

  mapped_mem = mmap(0, mmap_size, mmap_prots, mmap_flags, fd, 0);

  if (mapped_mem == MAP_FAILED)
  {
    std::cout << "Error mapping memory: "<< filename << std::endl;
    close(fd);
    return false;
  }

  ht.read_mmap(mapped_mem);

  close(fd);

  return true;
}

}

#endif

#endif
