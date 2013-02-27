/*  _________________________________________________________________________
 *
 *  MTGL: The MultiThreaded Graph Library
 *  Copyright (c) 2009 Sandia Corporation.
 *  This software is distributed under the BSD License.
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *  For more information, see the README file in the top MTGL directory.
 *  _________________________________________________________________________
 */

/****************************************************************************/
/*! \file test_sort.cpp

    \author Brad Mancke
    \author Greg Mackey (gemacke@sandia.gov)

    \date 6/4/2009

    Tests the sorting code.
*/
/****************************************************************************/

//#define DEBUG

#include <iostream>
#include <iomanip>

#include <mtgl/algorithm.hpp>
#include <mtgl/random.hpp>
#include <mtgl/util.hpp>
#include <mtgl/mtgl_io.hpp>

#ifdef USING_QTHREADS
#include <qthread.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace mtgl;

bool correctly_sorted(long* a, int size)
{
  int cnt = 0;

  for (int i = 1; i < size; i++)
  {
    if (a[i - 1] > a[i]) cnt++;
  }

  return cnt == 0;
}

bool is_number(char* c)
{
  int len = strlen(c);
  int i = 0;
  for ( ; i < len && isdigit(c[i]); ++i);

  return i == len;
}

int main(int argc, char* argv[])
{
#if defined _OPENMP || defined USING_QTHREADS
  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0] << " <filename | random_array_size> "
              << "<num threads>" << std::endl;
    exit(1);
  }
#else
  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " <filename | random_array_size>"
              << std::endl;
    exit(1);
  }
#endif

#ifdef USING_QTHREADS
  qthread_init(atoi(argv[2]));
#elif defined _OPENMP
  omp_set_num_threads(atoi(argv[2]));
#endif

  mt_timer mt;
  int size;
  long* arr;
  long* arr1;
  long* arr2;

  if (is_number(argv[1]))
  {
    #pragma mta fence
    mt.start();

    size = atoi(argv[1]);
    arr = (long*) malloc(size * sizeof(long));
    arr1 = (long*) malloc(size * sizeof(long));
    arr2 = (long*) malloc(size * sizeof(long));

    mt_lrand48(size, arr);

    int mod_size = size / 2;

    #pragma mta noalias *arr
    #pragma mta assert nodep
    for (int i = 0; i < size; i++)
    {
      arr[i] = arr[i] % mod_size + mod_size;
      arr1[i] = arr[i];
      arr2[i] = arr[i];
    }

    #pragma mta fence
    mt.stop();

    std::cout << "Array creation time: " << mt.getElapsedSeconds() << std::endl;
  }
  else
  {
    #pragma mta fence
    mt.start();

    arr = read_array<long>(argv[1], size); 
    arr1 = (long*) malloc(size * sizeof(long));
    arr2 = (long*) malloc(size * sizeof(long));

    #pragma mta noalias *arr
    #pragma mta assert nodep
    for (int i = 0; i < size; i++)
    {
      arr1[i] = arr[i];
      arr2[i] = arr[i];
    }

    #pragma mta fence
    mt.stop();

    std::cout << "Array read time: " << mt.getElapsedSeconds() << std::endl;
  }

  #pragma mta fence
  mt.start();
  sort(arr, size);
  #pragma mta fence
  mt.stop();

  std::cout << "Sorting one array time: " << mt.getElapsedSeconds()
            << std::endl;

  if (correctly_sorted(arr, size))
  {
    std::cout << "Array sorted correctly." << std::endl;
  }
  else
  {
    std::cout << "Array not sorted correctly." << std::endl;
  }

  #pragma mta fence
  mt.start();
  sort(arr1, size, arr2);
  #pragma mta fence
  mt.stop();

  std::cout << "Sorting two arrays time: " << mt.getElapsedSeconds()
            << std::endl;

  if (correctly_sorted(arr1, size))
  {
    std::cout << "Array1 sorted correctly." << std::endl;
  }
  else
  {
    std::cout << "Array1 not sorted correctly." << std::endl;
  }

  if (correctly_sorted(arr2, size))
  {
    std::cout << "Array2 sorted correctly." << std::endl;
  }
  else
  {
    std::cout << "Array2 not sorted correctly." << std::endl;
  }

  int cnt = 0;
  for (int i = 0; i < size; i++)
  {
    if (arr1[i] != arr2[i]) cnt++;
  }

  if (cnt)
  {
    std::cout << "Sorting two arrays didn't work out." << std::endl;
  }
  else
  {
    std::cout << "Sorting two arrays worked out." << std::endl;
  }

  free(arr);
  free(arr1);
  free(arr2);

  return 0;
}
