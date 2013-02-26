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
/*! \file random.hpp

    \brief This file contains random number generation functions that
           work on Unix, XMT, and Windows systems.  There are functions
           to return a single random number and functions to populate an
           array of random numbers.

           In most cases random() is the preferred function to generate a
           single random number.  It is thread-safe on the XMT, and on most
           systems it generates as good or better a random number as both
           rand() and the rand48() family of functions.  The one exception
           is Windows which does not have a random() function.  In this case
           we use rand().

           On the XMT the prand() family of functions is used to generate
           arrays of random numbers in parallel.  On all other systems, we 
           use a for loop and call the single random number generator to 
           populate the array.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 2/24/2011
*/
/****************************************************************************/

#ifndef MTGL_RANDOM_HPP
#define MTGL_RANDOM_HPP

#include <cstdlib>
#include <limits>
#include <iostream>

#include <mtgl/types.hpp>

#ifdef __MTA__
#include <mta_rng.h>
#endif

namespace mtgl {

inline
void mt_srand48(unsigned long seed)
{
#ifdef _WIN32
  srand(seed);
#else
  srand48(seed);
#endif
}

inline
long mt_lrand48()
{
#ifdef _WIN32
  return static_cast<long>(rand());
#else
  return lrand48();
#endif
}

inline
long mt_lrand48_64()
{
#if SIZEOF_LONG == 8
  return mt_lrand48() + (mt_lrand48() << 31);
#else
  return mt_lrand48();
#endif
}

inline
double mt_drand48()
{
#ifdef _WIN32
  return static_cast<double>(rand()) / RAND_MAX;
#else
  return drand48();
#endif
}

inline
void mt_lrand48(int size, long* vals)
{
#ifdef __MTA__
  prand_int(size, reinterpret_cast<int*>(vals));
#else
  for (int i = 0; i < size; ++i) vals[i] = mt_lrand48();
#endif
}

inline
void mt_lrand48_64(int size, long* vals)
{
#ifdef __MTA__
  prand_int(size, reinterpret_cast<int*>(vals));
#else
  for (int i = 0; i < size; ++i) vals[i] = mt_lrand48_64();
#endif
}

inline
void mt_drand48(int size, double* vals)
{
#ifdef __MTA__
  prand(size, vals);
#else
  for (int i = 0; i < size; ++i) vals[i] = mt_drand48();
#endif
}

inline
void mt_srandom(unsigned seed)
{
#ifdef _WIN32
  srand(seed);
#else
  srandom(seed);
#endif
}

inline
long mt_random()
{
#ifdef _WIN32
  return static_cast<long>(rand());
#else
  return random();
#endif
}

inline
long mt_random_64()
{
#if SIZEOF_LONG == 8
  return mt_random() + (mt_random() << 31);
#else
  return mt_random();
#endif
}

inline
double mt_drandom()
{
#ifdef _WIN32
  return static_cast<double>(rand()) / RAND_MAX;
#else
  return static_cast<double>(random()) / RAND_MAX;
#endif
}

inline
void mt_random(int size, long* vals)
{
#ifdef __MTA__
  prand_int(size, reinterpret_cast<int*>(vals));
#else
  for (int i = 0; i < size; ++i) vals[i] = mt_random();
#endif
}

inline
void mt_random_64(int size, long* vals)
{
#ifdef __MTA__
  prand_int(size, reinterpret_cast<int*>(vals));
#else
  for (int i = 0; i < size; ++i) vals[i] = mt_random_64();
#endif
}

inline
void mt_drandom(int size, double* vals)
{
#ifdef __MTA__
  prand(size, vals);
#else
  for (int i = 0; i < size; ++i) vals[i] = mt_drandom();
#endif
}

/// \brief Takes an array and permutes the elements
/// \param n The size of the array
/// \param array The array to be permuted.
template <typename T, typename size_type>
void random_permutation(size_type n, T* array)
{
  long* randVals = (long*) malloc(n * sizeof(long));
  mt_lrand48(n, randVals);

#ifdef __MTA__
  if (n > 10000)
  {
    #pragma mta assert parallel
    for (size_type i = 0; i < n; ++i)
    {
      size_type j = randVals[i] % n;
      if (i != j)
      {
        T x = mt_readfe(array[i]);
        T y = mt_readfe(array[j]);
        mt_write(array[i], y);
        mt_write(array[j], x);
      }
    }
  }
  else
  {
    for (size_type i = 0; i < n; ++i)
    {
      size_type j = randVals[i] % n;
      if (i != j)
      {
        T tmp = array[i];
        array[i] = array[j];
        array[j] = tmp;
      }
    }
  }
#else
  for (size_type i = 0; i < n; ++i)
  {
    size_type j = randVals[i] % ((unsigned int) n);
    if (i != j)
    {
      T tmp = array[i];
      array[i] = array[j];
      array[j] = tmp;
    }
  }
#endif

  free(randVals);
}

/// \brief Adapted from Kamesh Madduri's code for DIMACS shortest
//         path challenge.
template <typename T, typename size_type>
void random_permutation(size_type n, size_type m, T* srcs, T* dests)
{
  long* randVals = (long*) malloc(n * sizeof(long));
  mt_lrand48(n, randVals);

  T* perm = (T*) malloc(n * sizeof(T));

  for (size_type i = 0; i < n; ++i) perm[i] = i;

  #pragma mta assert parallel
  for (size_type i = 0; i < n; ++i)
  {
#ifdef __MTA__
    size_type j = randVals[i] % n;
#else
    size_type j = randVals[i] % ((unsigned int) n);
#endif

    if (i == j) continue;

    /* Swap perm[i] and perm[j] */
#ifdef __MTA__
    T x = mt_readfe(perm[i]);
    T y = mt_readfe(perm[j]);
    mt_write(perm[i], y);
    mt_write(perm[j], x);
#else
    T tmp = perm[i];
    perm[i] = perm[j];
    perm[j] = tmp;
#endif
  }

  #pragma mta assert nodep
  for (size_type i = 0; i < m; ++i)
  {
    srcs[i] = perm[srcs[i]];
    dests[i] = perm[dests[i]];
  }

  free(perm);
  free(randVals);
}

class lrand48_generator {
public:
#ifdef __MTA__
  inline
  lrand48_generator(int sz) : size(sz),
                              _store((long*) malloc(size * sizeof(long)))
  {
    prand_int(size, reinterpret_cast<int*>(_store));
  }

  inline ~lrand48_generator() { free(_store); }
  inline long operator[](int i) { return _store[i]; }

  inline void generate(int sz)
  {
    assert(sz <= size);
    prand_int(sz, reinterpret_cast<int*>(_store));
  }

private:
  int size;
  long* _store;
#else
  inline lrand48_generator(int sz) {}
  inline ~lrand48_generator() {}
  inline long operator[](int i) { return mt_lrand48(); }
  inline void generate(int sz) {}
#endif
};

class lrand48_64_generator {
public:
#ifdef __MTA__
  inline
  lrand48_64_generator(int sz) : size(sz),
                                 _store((long*) malloc(size * sizeof(long)))
  {
    prand_int(size, reinterpret_cast<int*>(_store));
  }

  inline ~lrand48_64_generator() { free(_store); }
  inline long operator[](int i) { return _store[i]; }

  inline void generate(int sz)
  {
    assert(sz <= size);
    prand_int(sz, reinterpret_cast<int*>(_store));
  }

private:
  int size;
  long* _store;
#else
  inline lrand48_64_generator(int sz) {}
  inline ~lrand48_64_generator() {}
  inline long operator[](int i) { return mt_lrand48_64(); }
  inline void generate(int sz) {}
#endif
};

class drand48_generator {
public:
#ifdef __MTA__
  inline
  drand48_generator(int sz) : size(sz),
                              _store((double*) malloc(size * sizeof(double)))
  {
    prand(size, _store);
  }

  inline ~drand48_generator() { free(_store); }
  inline double operator[](int i) { return _store[i]; }

  inline void generate(int sz)
  {
    assert(sz <= size);
    prand(sz, _store);
  }

private:
  int size;
  double* _store;
#else
  inline drand48_generator(int sz) {}
  inline ~drand48_generator() {}
  inline double operator[](int i) { return mt_drand48(); }
  inline void generate(int sz) {}
#endif
};

class random_generator {
public:
#ifdef __MTA__
  inline
  random_generator(int sz) : size(sz),
                             _store((long*) malloc(size * sizeof(long)))
  {
    prand_int(size, reinterpret_cast<int*>(_store));
  }

  inline ~random_generator() { free(_store); }
  inline long operator[](int i) { return _store[i]; }

  inline void generate(int sz)
  {
    assert(sz <= size);
    prand_int(sz, reinterpret_cast<int*>(_store));
  }

private:
  int size;
  long* _store;
#else
  inline random_generator(int sz) {}
  inline ~random_generator() {}
  inline long operator[](int i) { return mt_random(); }
  inline void generate(int sz) {}
#endif
};

class random_64_generator {
public:
#ifdef __MTA__
  inline
  random_64_generator(int sz) : size(sz),
                                _store((long*) malloc(size * sizeof(long)))
  {
    prand_int(size, reinterpret_cast<int*>(_store));
  }

  inline ~random_64_generator() { free(_store); }
  inline long operator[](int i) { return _store[i]; }

  inline void generate(int sz)
  {
    assert(sz <= size);
    prand_int(sz, reinterpret_cast<int*>(_store));
  }

private:
  int size;
  long* _store;
#else
  inline random_64_generator(int sz) {}
  inline ~random_64_generator() {}
  inline long operator[](int i) { return mt_random_64(); }
  inline void generate(int sz) {}
#endif
};

class drandom_generator {
public:
#ifdef __MTA__
  inline
  drandom_generator(int sz) : size(sz),
                              _store((double*) malloc(size * sizeof(double)))
  {
    prand(size, _store);
  }

  inline ~drandom_generator() { free(_store); }
  inline double operator[](int i) { return _store[i]; }

  inline void generate(int sz)
  {
    assert(sz <= size);
    prand(sz, _store);
  }

private:
  int size;
  double* _store;
#else
  inline drandom_generator(int sz) {}
  inline ~drandom_generator() {}
  inline double operator[](int i) { return mt_drandom(); }
  inline void generate(int sz) {}
#endif
};

}

#endif
