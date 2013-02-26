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
/*! \file numeric.hpp

    Mimics some of the functionality of numeric in the standard template
    library.

    \author Eric Goodman (elgoodm@sandia.gov)
    \author Greg Mackey (gemacke@sandia.gov)

    \date 11/2011
*/
/****************************************************************************/

#ifndef MTGL_NUMERIC_HPP
#define MTGL_NUMERIC_HPP

#include <mtgl/partitioning.hpp>
#include <mtgl/util.hpp>

#ifdef USING_QTHREADS
#include <qthread/qloop.hpp>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace mtgl {

#ifdef USING_QTHREADS
namespace detail {

template <typename InputIterator, typename T>
class accum_loop {
public:
  accum_loop(const InputIterator& f, T& t) : first(f), total(t) {}

  void operator()(size_t start, size_t stop)
  {
    T my_total = 0;

    InputIterator my_first = first + start;
    InputIterator my_last = first + stop;

    for ( ; my_first != my_last; ++my_first) my_total += *my_first;

    mt_incr(total, my_total);
  }

private:
  const InputIterator& first;
  T& total;
};

}
#endif

template <typename T, typename InputIterator>
T accumulate(const InputIterator& first, const InputIterator& last, T init)
{
  T total = init;
  unsigned long size = last - first;

#ifdef USING_QTHREADS
  detail::accum_loop<InputIterator, T> accum_obj(first, total);
  qt_loop_balance(0, size, accum_obj);
#else
  unsigned long thread_id = 0;
  unsigned long num_threads = 1;

  #pragma mta for all streams thread_id of num_threads
  #pragma mta use 100 streams
#ifdef _OPENMP
  #pragma omp parallel default(shared)
#endif
  {
#ifdef _OPENMP
    unsigned long thread_id = omp_get_thread_num();
    unsigned long num_threads = omp_get_num_threads();
#endif

    unsigned long start = begin_block_range(size, thread_id, num_threads);
    unsigned long end = end_block_range(size, thread_id, num_threads);

    T my_total = 0;

    InputIterator my_first = first + start;
    InputIterator my_last = first + end;

    for ( ; my_first != my_last; ++my_first) my_total += *my_first;

    mt_incr(total, my_total);
  }
#endif

  #pragma mta fence
  return total;
}

}

#endif
