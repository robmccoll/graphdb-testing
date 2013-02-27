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
/*! \file test_random.cpp

    \brief Test the random functions in random.hpp.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 2/2011
*/
/****************************************************************************/

#include <iostream>
#include <limits>
#include <iomanip>

#include <mtgl/random.hpp>

using namespace mtgl;

int main()
{
  const int num_vals = 70;

  std::cout << "          RAND_MAX: " << RAND_MAX << std::endl << std::endl;

  std::cout << "random functions:" << std::endl;
  long v_long  = mt_random();
  std::cout << "    long: " << v_long << std::endl;

  long v_long_64  = mt_random_64();
  std::cout << "big long: " << v_long_64 << std::endl;

  double v_double  = mt_drandom();
  std::cout << "  double: "
            << std::setprecision(std::numeric_limits<double>::digits)
            << v_double << std::endl << std::endl;

  std::cout << "rand48 functions:" << std::endl;
  v_long  = mt_lrand48();
  std::cout << "    long: " << v_long << std::endl;

  v_long_64  = mt_lrand48_64();
  std::cout << "big long: " << v_long_64 << std::endl;

  v_double  = mt_drand48();
  std::cout << "  double: "
            << std::setprecision(std::numeric_limits<double>::digits)
            << v_double << std::endl << std::endl;

  long* lvals = (long*) malloc(sizeof(long) * num_vals);
  mt_lrand48(num_vals, lvals);

  std::cout << "lrand48:" << std::endl;
  for (int i = 0; i < num_vals; ++i)  std::cout << lvals[i] << std::endl;
  std::cout << std::endl;

  free(lvals);

  double* dvals = (double*) malloc(sizeof(double) * num_vals);
  mt_drand48(num_vals, dvals);

  std::cout << "drand48:" << std::endl;
  for (int i = 0; i < num_vals; ++i)
  {
    std::cout << std::setprecision(std::numeric_limits<double>::digits)
              << dvals[i] << std::endl;
  }

  free(dvals);
  return 0;
}
