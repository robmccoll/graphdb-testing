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
/*! \file mtgl_subiso_triangles2.cpp

    \brief Defines the entry point for the console application.

    \author Greg Mackey (gemacke@sandia.gov)
    \author Jon Berry (jberry@sandia.gov)

    \date 2/9/2011
*/
/****************************************************************************/

#include "stdafx.h"

#define _CRTDBG_MAP_ALLOC

#include <cstdlib>
#include <crtdbg.h>
#include <test/subiso_triangles.cpp>

int _tmain(int argc, TCHAR* argv[])
{
  _CrtSetDbgFlag (_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);

  char** argv2 = (char**) malloc(sizeof(char*) * argc);
  for (int i = 0; i < argc; ++i)
  {
    argv2[i] = (char*) malloc(sizeof(char) * 256);
    wcstombs(argv2[i], argv[i], 256);
  }

  printf("n: %d, p: %f, s: %d\n", atoi(argv2[1]), atof(argv2[2]),
         atoi(argv2[3]));

  subiso_main(argc, argv2);

  for (int i = 0; i < argc; ++i) free(argv2[i]);
  free(argv2);

  getchar();
}
