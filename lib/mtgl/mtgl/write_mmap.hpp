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
/*! \file write_mmap.hpp

    \brief Writes a graph to an mmapped memory region.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 11/11/2010
*/
/****************************************************************************/

#ifndef MTGL_WRITE_MMAP_HPP
#define MTGL_WRITE_MMAP_HPP

#ifndef _WIN32

#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>

#include <mtgl/mtgl_adapter.hpp>

namespace mtgl {

template <typename Graph>
bool write_mmap(Graph& g, const char* filename)
{
  typedef typename Graph::size_type size_type;

  int mmap_prots = PROT_READ | PROT_WRITE;
  int mmap_flags = MAP_SHARED;

#ifdef __MTA__
  mmap_flags |= MAP_ANON;
  int fd = open(filename, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
#else
  int fd = shm_open(filename, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
#endif

  if (fd < 0)
  {
    std::cout << "Error opening file: "<< filename << std::endl;
    return;
  }

  size_type mapped_size = g.get_mmap_size();

#ifndef __MTA__
  if (ftruncate(fd, mapped_size) != 0)
  {
    std::cout << "Error resizing memory region: "<< filename << "    size: "
              << mapped_size << std::endl;
    return;
  }
#endif

  void* mapped_mem = mmap(0, mapped_size, mmap_prots, mmap_flags, fd, 0);

  if (mapped_mem == MAP_FAILED)
  {
    close(fd);

#ifdef __MTA__
    remove(filename);
#else
    shm_unlink(filename);
#endif

    std::cout << "Error mapping memory: "<< filename << std::endl;
    return false;
  }

  g.write_mmap(mapped_mem);

  close(fd);

  return true;
}

}

#endif 

#endif
