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
/*! \file test_shared_ptr.cpp

    \brief Tests the functionality of the shared_ptr.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 12/8/2010
*/
/****************************************************************************/

#include <cstdio>

#include <mtgl/shared_ptr.hpp>

using namespace mtgl;

int main()
{
  shared_ptr<int> p(new int());
  shared_ptr<int> p2;

  if (p) printf("p: %d\n", *p);
  if (p2) printf("p2: %d\n", *p2);

  *p = 2;

  if (p) printf("p: %d\n", *p);
  if (p2) printf("p2: %d\n", *p2);

  p2 = p;

  if (p) printf("p: %d\n", *p);
  if (p2) printf("p2: %d\n", *p2);

  p2.reset();

  if (p) printf("p: %d\n", *p);
  if (p2) printf("p2: %d\n", *p2);

  return 0;
}
