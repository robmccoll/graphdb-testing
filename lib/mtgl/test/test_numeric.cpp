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
/*! \file test_numeric.hpp

    Tests functionality in numeric.hpp

    \author Eric Goodman (elgoodm@sandia.gov)

    \date 11/2011

*/
/****************************************************************************/


#include <mtgl/numeric.hpp>
#include <mtgl/util.hpp>

#ifdef USING_QTHREADS
#include <qthread.h>
#endif

using namespace mtgl;

int main(int argc, const char* argv[])
{
  typedef int size_type;

  size_type cur_arg_index = 1;  ///>Used to process the command line arguments
  mt_timer timer;               ///>For performance info
  size_type n = 1000;

  const char* help_string =
    "This tests functionality in numeric.hpp.\n"
    "The format of the command is the following:\n\n"
    "test_numeric\n"
#ifdef USING_QTHREADS
    "\t-nt <number>    : Number of qthreads to use.\n"
#endif
    "\t-size <number>  : The size of the tests.\n";

  if (argc < 3)
  {
    std::cout << help_string << std::endl;
    return 0;
  }

#ifdef USING_QTHREADS
  unsigned long num_threads = 0;
#endif

  while (cur_arg_index < argc)
  {
    if (strcmp(argv[cur_arg_index], "-size") == 0)
    {
      ++cur_arg_index;
      n = strtoul(argv[cur_arg_index], NULL, 10);
    }
    else if (strcmp(argv[cur_arg_index], "--help") == 0 ||
             strcmp(argv[cur_arg_index], "-help") == 0)
    {
      std::cout << help_string << std::endl;
      return 0;
    }
#ifdef USING_QTHREADS
    else if (strcmp(argv[cur_arg_index], "-nt") == 0)
    {
      ++cur_arg_index;
      num_threads = strtoul(argv[cur_arg_index], NULL, 10);
    }
#endif
    else
    {
      std::cout << "Unrecognized option: " << argv[cur_arg_index] << std::endl;
    }

    ++cur_arg_index;
  }

#ifdef USING_QTHREADS
  qthread_init(num_threads);
#endif

  std::cout << "Testing numeric.hpp." << std::endl
            << "Testing accumulate()." << std::endl;

  timer.start();

  size_type* array = (size_type*) malloc(sizeof(size_type) * n);

#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for (size_type i = 0; i < n; ++i) array[i] = 1;

  #pragma mta fence
  timer.stop();
  std::cout << "\tArray setup time: " << timer.getElapsedSeconds()
            << std::endl;

  timer.start();

  size_type answer = accumulate(array, array + n, static_cast<size_type>(0));

  #pragma mta fence
  timer.stop();
  std::cout << "\tAccumulate time: " << timer.getElapsedSeconds()
            << std::endl;

  std::cout << "\tNum entries: " << n << std::endl;

  if (answer != n)
  {
    std::cout << "\tERROR in accumulate.  Expected " << n << ", got " << answer
              << std::endl;
  }
  else
  {
    std::cout << "\tResult: " << answer << std::endl;
  }

  free(array);
}
