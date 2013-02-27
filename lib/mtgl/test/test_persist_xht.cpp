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
/*! \file test_persist_xht.cpp

    \brief Test to see if mmapping an xmt_hash_table works.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 5/25/2011
*/
/****************************************************************************/

#include <string>
#include <cstdio>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>

#include <mtgl/xmt_hash_table.hpp>

using namespace mtgl;

int main(int argc, char* argv[])
{
  xmt_hash_table<int, int> ht(10);

  for (int i = 0; i < 10; ++i)
  {
    ht.insert(i, i + 10);
  }
  ht.print();

  char cur_path[FILENAME_MAX];
  getcwd(cur_path, sizeof(cur_path));

  char file_name[] = "xht.mmap";

  int mmap_prots = PROT_READ | PROT_WRITE;
  int mmap_flags = MAP_SHARED;

#ifdef __MTA__
  mmap_flags |= MAP_ANON;
  int fd = open(file_name, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
#else
  int fd = shm_open(file_name, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
#endif

  if (fd < 0)
  {
    std::cout << "Error opening file: "<< file_name << std::endl;
    return 0;
  }

  unsigned long mapped_size = ht.get_mmap_size();

#ifndef __MTA__
  if (ftruncate(fd, mapped_size) != 0)
  {
    std::cout << "Error resizing memory region: "<< file_name << std::endl;
    return 0;
  }
#endif

  void* mapped_mem = mmap(0, mapped_size, mmap_prots, mmap_flags, fd, 0);

  if (mapped_mem == MAP_FAILED)
  {
    close(fd);

#ifdef __MTA__
    remove(file_name);
#else
    shm_unlink(file_name);
#endif

    std::cout << "Error mapping memory: "<< file_name << std::endl;
    return 0;
  }

  ht.write_mmap(mapped_mem);

  close(fd);

  std::cout << std::endl << "<Press enter to exit.>";
  std::string input;
  std::getline(std::cin, input);

#ifdef __MTA__
    remove(file_name);
#else
    shm_unlink(file_name);
#endif

  return 0;
}
