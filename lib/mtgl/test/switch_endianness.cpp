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
/*! \file switch_endianness.cpp

    \brief Reads in an array of unsigned longs and writes it back out again
           with the endianness of each of the array entries reversed.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 9/14/2011
*/
/****************************************************************************/

#define DEBUG

#include <mtgl/mtgl_io.hpp>
#include <mtgl/util.hpp>

using namespace mtgl;

int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    fprintf(stderr, "Usage: %s <filename>\n", argv[0]);

    exit(1);
  }

  char* filename = argv[1];

  unsigned long array_size;
  unsigned long* array = read_array<unsigned long>(filename, array_size);

  byte_swap(array, array_size, sizeof(unsigned long));

  write_array(filename, array, array_size);

  free(array);

  return 0;
}
