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
/*! \file algorithm.hpp

    Mimics some of the functionality of algorithm in the standard template
    library.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 01/2012
*/
/****************************************************************************/

#ifndef MTGL_ALGORITHM_HPP
#define MTGL_ALGORITHM_HPP

#include <algorithm>
#include <iostream>
#include <iomanip>

#include <mtgl/merge_sort.hpp>
#include <mtgl/partitioning.hpp>

#ifdef USING_QTHREADS
#include <qthread/qloop.hpp>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace mtgl {

namespace detail {

#if defined _OPENMP || defined USING_QTHREADS
template <typename T, typename T2>
class sort_two_array_iterator {
public:
  struct pair_ref;

  struct pair_val {
    pair_val(T f, T2 s) : first(f), second(s) {}

    pair_val(const pair_ref& rhs) : first(rhs.first), second(rhs.second) {}

    pair_val& operator=(const pair_val& rhs)
    {
      first = rhs.first;
      second = rhs.second;
      return *this;
    }

    pair_val& operator=(const pair_ref& rhs)
    {
      first = rhs.first;
      second = rhs.second;
      return *this;
    }

    T first;
    T2 second;
  };

  struct pair_ref {
    pair_ref(T& f, T2& s) : first(f), second(s) {}

    pair_ref(const pair_val& rhs) : first(rhs.first), second(rhs.second) {}

    pair_ref& operator=(const pair_ref& rhs)
    {
      first = rhs.first;
      second = rhs.second;
      return *this;
    }

    pair_ref& operator=(const pair_val& rhs)
    {
      first = rhs.first;
      second = rhs.second;
      return *this;
    }

    T& first;
    T2& second;
  };

  typedef ptrdiff_t difference_type;
  typedef pair_val value_type;
  typedef pair_val* pointer;
  typedef pair_ref reference;
  typedef std::random_access_iterator_tag iterator_category;

  sort_two_array_iterator(T* aa, T2* bb, long pp) : a(aa), b(bb), p(pp) {}

  sort_two_array_iterator(const sort_two_array_iterator& rhs) :
    a(rhs.a), b(rhs.b), p(rhs.p) {}

  sort_two_array_iterator& operator++()
  {
    ++p;
    return *this;
  }

  sort_two_array_iterator operator++(int)
  {
    sort_two_array_iterator temp(*this);
    ++p;
    return temp;
  }

  sort_two_array_iterator& operator--()
  {
    --p;
    return *this;
  }

  sort_two_array_iterator operator--(int)
  {
    sort_two_array_iterator temp(*this);
    --p;
    return temp;
  }

  sort_two_array_iterator operator+(ptrdiff_t n)
  {
    sort_two_array_iterator temp(*this);
    temp.p += n;
    return temp;
  }

  sort_two_array_iterator operator-(ptrdiff_t n)
  {
    sort_two_array_iterator temp(*this);
    temp.p -= n;
    return temp;
  }

  ptrdiff_t operator-(sort_two_array_iterator& i)
  {
    return p - i.p;
  }

  reference operator*() { return reference(a[p], b[p]); }

  bool operator==(const sort_two_array_iterator& rhs) const
  { return p == rhs.p; }
  bool operator!=(const sort_two_array_iterator& rhs) const
  { return p != rhs.p; }
  bool operator<(const sort_two_array_iterator& rhs) const
  { return p < rhs.p; }
  bool operator>(const sort_two_array_iterator& rhs) const
  { return p > rhs.p; }
  bool operator<=(const sort_two_array_iterator& rhs) const
  { return p <= rhs.p; }
  bool operator>=(const sort_two_array_iterator& rhs) const
  { return p >= rhs.p; }

private:
  T* a;
  T2* b;
  long p;
};

template <typename T, typename T2, typename Comparator>
class sort_two_array_iterator_comp {
public:
  typedef typename sort_two_array_iterator<T, T2>::value_type VT;
  typedef typename sort_two_array_iterator<T, T2>::reference RT;

  sort_two_array_iterator_comp(Comparator c) : comp(c) {}

  bool
  operator()(const VT& t1, const VT& t2) { return comp(t1.first, t2.first); }

  bool
  operator()(const RT& t1, const VT& t2) { return comp(t1.first, t2.first); }

  bool
  operator()(const VT& t1, const RT& t2) { return comp(t1.first, t2.first); }

  bool
  operator()(const RT& t1, const RT& t2) { return comp(t1.first, t2.first); }

private:
  Comparator comp;
};

#ifdef USING_QTHREADS
template <typename T, typename Comparator>
class thread_sort_qt {
public:
  thread_sort_qt(T* ar, long* p, Comparator c) :
    a(ar), partitions(p), comp(c) {}

  void operator()(const size_t start, const size_t stop)
  {
    long start_pos = partitions[start];
    long end_pos = partitions[start + 1];

    std::sort(a + start_pos, a + end_pos, comp);
  }

private:
  T* a;
  long* partitions;
  Comparator comp;
};

template <typename T, typename Comparator>
class thread_merge_qt {
public:
  thread_merge_qt(T* ar, T* br, long* p, Comparator c) :
    a(ar), b(br), partitions(p), comp(c) {}

  void operator()(const size_t start, const size_t stop)
  {
    if (start != stop)
    {
      long left_pos = partitions[start * 2];
      long left_end = partitions[start * 2 + 1];
      long right_pos = partitions[start * 2 + 1];
      long right_end = partitions[(start + 1) * 2];
      long b_pos = partitions[start * 2];

      // Merge the left and right halves until one of them is empty.
      while (left_pos < left_end && right_pos < right_end)
      {
        if (comp(a[left_pos], a[right_pos]))
        {
          b[b_pos++] = a[left_pos++];
        }
        else
        {
          b[b_pos++] = a[right_pos++];
        }
      }

      // There will be a tail left in either the left or right half.  Copy
      // the tail to b.
      if (left_pos < left_end)
      {
        // Copy the tail in the left half to b.
        for (long j = left_pos; j < left_end; ++j)
        {
          b[b_pos + (j - left_pos)] = a[j];
        }
      }
      else
      {
        // Copy the tail in the right half to b.
        for (long j = right_pos; j < right_end; ++j)
        {
          b[b_pos + (j - right_pos)] = a[j];
        }
      }
    }
  }

private:
  T* a;
  T* b;
  long* partitions;
  Comparator comp;
};

template <typename T, typename T2, typename Comparator>
class thread_sort_qt2 {
public:
  thread_sort_qt2(T* ar, T2* ar2, long* p, Comparator c) :
    a(ar), a2(ar2), partitions(p), comp(c) {}

  void operator()(const size_t start, const size_t stop)
  {
    long start_pos = partitions[start];
    long end_pos = partitions[start + 1];

    sort_two_array_iterator<T, T2> start_iter(a, a2, start_pos);
    sort_two_array_iterator<T, T2> end_iter(a, a2, end_pos);

    std::sort(start_iter, end_iter,
              sort_two_array_iterator_comp<T, T2, Comparator>(comp));
  }

private:
  T* a;
  T2* a2;
  long* partitions;
  Comparator comp;
};

template <typename T, typename T2, typename Comparator>
class thread_merge_qt2 {
public:
  thread_merge_qt2(T* ar, T* br, T2* ar2, T2* br2, long* p, Comparator c) :
    a(ar), b(br), a2(ar2), b2(br2), partitions(p), comp(c) {}

  void operator()(const size_t start, const size_t stop)
  {
    if (start != stop)
    {
      long left_pos = partitions[start * 2];
      long left_end = partitions[start * 2 + 1];
      long right_pos = partitions[start * 2 + 1];
      long right_end = partitions[(start + 1) * 2];
      long b_pos = partitions[start * 2];

      // Merge the left and right halves until one of them is empty.
      while (left_pos < left_end && right_pos < right_end)
      {
        if (comp(a[left_pos], a[right_pos]))
        {
          b[b_pos] = a[left_pos];
          b2[b_pos++] = a2[left_pos++];
        }
        else
        {
          b[b_pos] = a[right_pos];
          b2[b_pos++] = a2[right_pos++];
        }
      }

      // There will be a tail left in either the left or right half.  Copy
      // the tail to b.
      if (left_pos < left_end)
      {
        // Copy the tail in the left half to b.
        for (long j = left_pos; j < left_end; ++j)
        {
          b[b_pos + (j - left_pos)] = a[j];
          b2[b_pos + (j - left_pos)] = a2[j];
        }
      }
      else
      {
        // Copy the tail in the right half to b.
        for (long j = right_pos; j < right_end; ++j)
        {
          b[b_pos + (j - right_pos)] = a[j];
          b2[b_pos + (j - right_pos)] = a2[j];
        }
      }
    }
  }

private:
  T* a;
  T* b;
  T2* a2;
  T2* b2;
  long* partitions;
  Comparator comp;
};

#ifndef MTGL_ARR_COPY_LOOP
#define MTGL_ARR_COPY_LOOP
template <typename T>
class arr_copy_loop {
public:
  arr_copy_loop(T* d, T* s) : dest(d), src(s) {}

  void operator()(const size_t start, const size_t stop)
  {
    for (size_t i = start; i != stop; ++i) dest[i] = src[i];
  }

private:
  T* dest;
  T* src;
};
#endif
#endif

template <typename T, typename Comparator>
void omp_qt_sort(T* array, long size, Comparator comp)
{
  // a starts out pointing to the original array, and b starts out pointing to
  // the temporary memory.  a and b are going to be swapped back and forth by
  // the algorithm.  At the beginning of each iteration, a holds the unmerged
  // arrays.  The merged results are put into b.  Then the pointers for a and
  // b are swapped, so the result is pointed to by a at the end.
  T* new_array = (T*) malloc(size * sizeof(T));
  T* a = array;
  T* b = new_array;

#ifdef DEBUG
  std::cout << "Original:";
  for (int i = 0; i < size; i++)
  {
    std::cout << std::fixed << std::setw(5) << array[i];
  }
  std::cout << std::endl;
#endif

#ifdef USING_QTHREADS
  long num_threads = qthread_num_shepherds();
#else
  long num_threads = omp_get_max_threads();
#endif
  long* partitions = (long*) malloc((num_threads + 1) * sizeof(long));

  // Determine the partitions for each thread.
  for (long i = 0; i < num_threads; ++i)
  {
    partitions[i] = begin_block_range(size, i, num_threads);
  }
  partitions[num_threads] = end_block_range(size, num_threads - 1, num_threads);

  // Each thread sorts an equal sized block of the data.
#ifdef USING_QTHREADS
  detail::thread_sort_qt<T, Comparator> tsq(a, partitions, comp);
  qt_loop_balance(0, num_threads, tsq);
#else
  #pragma omp parallel
  {
    long thread_id = omp_get_thread_num();
    long start_pos = partitions[thread_id];
    long end_pos = partitions[thread_id + 1];

    std::sort(a + start_pos, a + end_pos, comp);
  }
#endif

#ifdef DEBUG
  std::cout << "  Pieces:";
  for (long i = 0; i < size; i++)
  {
    std::cout << std::fixed << std::setw(5) << a[i];
  }
  std::cout << std::endl;
#endif

  // Now, we need to merge the sorted lists.  If there are N threads, we use
  // N/2 threads to merge the N sorted lists.  Then we use N/4 threads to
  // merge the N/2 sorted lists and so on until we have a single sorted list.
  while (num_threads > 1)
  {
    num_threads /= 2;

#ifdef USING_QTHREADS
  detail::thread_merge_qt<T, Comparator> tmq(a, b, partitions, comp);
  qt_loop_balance(0, num_threads, tmq);
#else
    #pragma omp parallel
    {
      long thread_id = omp_get_thread_num();

      if (thread_id < num_threads)
      {
        long left_pos = partitions[thread_id * 2];
        long left_end = partitions[thread_id * 2 + 1];
        long right_pos = partitions[thread_id * 2 + 1];
        long right_end = partitions[(thread_id + 1) * 2];
        long b_pos = partitions[thread_id * 2];

        // Merge the left and right halves until one of them is empty.
        while (left_pos < left_end && right_pos < right_end)
        {
          if (comp(a[left_pos], a[right_pos]))
          {
            b[b_pos++] = a[left_pos++];
          }
          else
          {
            b[b_pos++] = a[right_pos++];
          }
        }

        // There will be a tail left in either the left or right half.  Copy
        // the tail to b.
        if (left_pos < left_end)
        {
          // Copy the tail in the left half to b.
          for (long j = left_pos; j < left_end; ++j)
          {
            b[b_pos + (j - left_pos)] = a[j];
          }
        }
        else
        {
          // Copy the tail in the right half to b.
          for (long j = right_pos; j < right_end; ++j)
          {
            b[b_pos + (j - right_pos)] = a[j];
          }
        }
      }
    }
#endif

#ifdef DEBUG
  std::cout << "  Merged:";
  for (long i = 0; i < size; i++)
  {
    std::cout << std::fixed << std::setw(5) << b[i];
  }
  std::cout << std::endl;
#endif

    for (long i = 1; i <= num_threads; ++i)
    {
      partitions[i] = partitions[i * 2];
    }

    // Swap a and b.
    T* tmp = a;
    a = b;
    b = tmp;
  }

  // The answer is pointed to by a.  If a ends up pointing to the temporary
  // memory, copy the result to the passed in array.
  if (a != array)
  {
#ifdef USING_QTHREADS
    detail::arr_copy_loop<T> a_cpy(array, a);
    qt_loop_balance(0, size, a_cpy);
#else
    #pragma omp parallel for
    for (long i = 0; i < size; ++i) array[i] = a[i];
#endif
  }

  free(partitions);
  free(new_array);

#ifdef DEBUG
  std::cout << "  Sorted:";
  for (long i = 0; i < size; i++)
  {
    std::cout << std::fixed << std::setw(5) << array[i];
  }
  std::cout << std::endl;
#endif
}

template <typename T, typename T2, typename Comparator>
void omp_qt_sort(T* array, long size, T2* array2, Comparator comp)
{
  // a starts out pointing to the original array, and b starts out pointing to
  // the temporary memory.  a and b are going to be swapped back and forth by
  // the algorithm.  At the beginning of each iteration, a holds the unmerged
  // arrays.  The merged results are put into b.  Then the pointers for a and
  // b are swapped, so the result is pointed to by a at the end.
  T* new_array = (T*) malloc(size * sizeof(T));
  T* a = array;
  T* b = new_array;

  T2* new_array2 = (T2*) malloc(size * sizeof(T2));
  T2* a2 = array2;
  T2* b2 = new_array2;

#ifdef DEBUG
  std::cout << "Original a:";
  for (int i = 0; i < size; i++)
  {
    std::cout << std::fixed << std::setw(5) << array[i];
  }
  std::cout << std::endl;
  std::cout << "Original b:";
  for (int i = 0; i < size; i++)
  {
    std::cout << std::fixed << std::setw(5) << array2[i];
  }
  std::cout << std::endl;
#endif

#ifdef USING_QTHREADS
  long num_threads = qthread_num_shepherds();
#else
  long num_threads = omp_get_max_threads();
#endif
  long* partitions = (long*) malloc((num_threads + 1) * sizeof(long));

  // Determine the partitions for each thread.
  for (long i = 0; i < num_threads; ++i)
  {
    partitions[i] = begin_block_range(size, i, num_threads);
  }
  partitions[num_threads] = end_block_range(size, num_threads - 1, num_threads);

  // Each thread sorts an equal sized block of the data.
#ifdef USING_QTHREADS
  detail::thread_sort_qt2<T, T2, Comparator> tsq2(a, a2, partitions, comp);
  qt_loop_balance(0, num_threads, tsq2);
#else
  #pragma omp parallel
  {
    long thread_id = omp_get_thread_num();
    long start_pos = partitions[thread_id];
    long end_pos = partitions[thread_id + 1];

    sort_two_array_iterator<T, T2> start_iter(a, a2, start_pos);
    sort_two_array_iterator<T, T2> end_iter(a, a2, end_pos);

    std::sort(start_iter, end_iter,
              sort_two_array_iterator_comp<T, T2, Comparator>(comp));
  }
#endif

#ifdef DEBUG
  std::cout << "  Pieces a:";
  for (long i = 0; i < size; i++)
  {
    std::cout << std::fixed << std::setw(5) << a[i];
  }
  std::cout << std::endl;
  std::cout << "  Pieces b:";
  for (long i = 0; i < size; i++)
  {
    std::cout << std::fixed << std::setw(5) << a2[i];
  }
  std::cout << std::endl;
#endif

  // Now, we need to merge the sorted lists.  If there are N threads, we use
  // N/2 threads to merge the N sorted lists.  Then we use N/4 threads to
  // merge the N/2 sorted lists and so on until we have a single sorted list.
  while (num_threads > 1)
  {
    num_threads /= 2;

#ifdef USING_QTHREADS
  detail::thread_merge_qt2<T, T2, Comparator>
    tmq2(a, b, a2, b2, partitions, comp);
  qt_loop_balance(0, num_threads, tmq2);
#else
    #pragma omp parallel
    {
      long thread_id = omp_get_thread_num();

      if (thread_id < num_threads)
      {
        long left_pos = partitions[thread_id * 2];
        long left_end = partitions[thread_id * 2 + 1];
        long right_pos = partitions[thread_id * 2 + 1];
        long right_end = partitions[(thread_id + 1) * 2];
        long b_pos = partitions[thread_id * 2];

        // Merge the left and right halves until one of them is empty.
        while (left_pos < left_end && right_pos < right_end)
        {
          if (comp(a[left_pos], a[right_pos]))
          {
            b[b_pos] = a[left_pos];
            b2[b_pos++] = a2[left_pos++];
          }
          else
          {
            b[b_pos] = a[right_pos];
            b2[b_pos++] = a2[right_pos++];
          }
        }

        // There will be a tail left in either the left or right half.  Copy
        // the tail to b.
        if (left_pos < left_end)
        {
          // Copy the tail in the left half to b.
          for (long j = left_pos; j < left_end; ++j)
          {
            b[b_pos + (j - left_pos)] = a[j];
            b2[b_pos + (j - left_pos)] = a2[j];
          }
        }
        else
        {
          // Copy the tail in the right half to b.
          for (long j = right_pos; j < right_end; ++j)
          {
            b[b_pos + (j - right_pos)] = a[j];
            b2[b_pos + (j - right_pos)] = a2[j];
          }
        }
      }
    }
#endif

#ifdef DEBUG
    std::cout << "  Merged a:";
    for (long i = 0; i < size; i++)
    {
      std::cout << std::fixed << std::setw(5) << b[i];
    }
    std::cout << std::endl;
    std::cout << "  Merged b:";
    for (long i = 0; i < size; i++)
    {
      std::cout << std::fixed << std::setw(5) << b2[i];
    }
    std::cout << std::endl;
#endif

    for (long i = 1; i <= num_threads; ++i)
    {
      partitions[i] = partitions[i * 2];
    }

    // Swap a and b.
    T* tmp = a;
    a = b;
    b = tmp;

    // Swap a2 and b2.
    T2* tmp2 = a2;
    a2 = b2;
    b2 = tmp2;
  }

  // The answer is pointed to by a.  If a ends up pointing to the temporary
  // memory, copy the result to the passed in array.
  if (a != array)
  {
#ifdef USING_QTHREADS
    detail::arr_copy_loop<T> a_cpy(array, a);
    qt_loop_balance(0, size, a_cpy);
    detail::arr_copy_loop<T2> a2_cpy(array2, a2);
    qt_loop_balance(0, size, a2_cpy);
#else
    #pragma omp parallel for
    for (long i = 0; i < size; ++i) array[i] = a[i];

    #pragma omp parallel for
    for (long i = 0; i < size; ++i) array2[i] = a2[i];
#endif
  }

  free(partitions);
  free(new_array);
  free(new_array2);

#ifdef DEBUG
  std::cout << "  Sorted a:";
  for (long i = 0; i < size; i++)
  {
    std::cout << std::fixed << std::setw(5) << array[i];
  }
  std::cout << std::endl;
  std::cout << "  Sorted b:";
  for (long i = 0; i < size; i++)
  {
    std::cout << std::fixed << std::setw(5) << array2[i];
  }
  std::cout << std::endl;
#endif
}
#endif
}

template <typename T>
inline
void sort(T* array, long size)
{
#if defined _OPENMP || defined USING_QTHREADS
  detail::omp_qt_sort(array, size, std::less<T>());
#else
  merge_sort(array, size, std::less<T>());
#endif
}

template <typename T, typename Comparator>
inline
void sort(T* array, long size, Comparator comp)
{
#if defined _OPENMP || defined USING_QTHREADS
  detail::omp_qt_sort(array, size, comp);
#else
  merge_sort(array, size, comp);
#endif
}

template <typename T, typename T2>
inline
void sort(T* array, long size, T2* array2)
{
#if defined _OPENMP || defined USING_QTHREADS
  detail::omp_qt_sort(array, size, array2, std::less<T>());
#else
  merge_sort(array, size, array2, std::less<T>());
#endif
}

template <typename T, typename T2, typename Comparator>
inline
void sort(T* array, long size, T2* array2, Comparator comp)
{
#if defined _OPENMP || defined USING_QTHREADS
  detail::omp_qt_sort(array, size, array2, comp);
#else
  merge_sort(array, size, array2, comp);
#endif
}

template <typename T, typename Comparator>
void insertion_sort(T* array, long size, Comparator comp)
{
  for (long i = 1; i < size; i++)
  {
    T value = array[i];
    long j = i;

    if (comp(array[i], array[0]))
    {
      for ( ; j > 0; --j) array[j] = array[j - 1];

      array[0] = value;
    }
    else
    {
      for ( ; comp(value, array[j - 1]); --j) array[j] = array[j - 1];

      array[j] = value;
    }
  }
}

template <typename T>
void insertion_sort(T* array, long size)
{
  insertion_sort(array, size, std::less<T>());
}

template <typename T, typename T2, typename Comparator>
void insertion_sort(T* array, long size, T2* array2, Comparator comp)
{
  for (long i = 1; i < size; i++)
  {
    T value = array[i];
    T2 value2 = array2[i];
    long j = i;

    if (comp(array[i], array[0]))
    {
      for ( ; j > 0; --j)
      {
        array[j] = array[j - 1];
        array2[j] = array2[j - 1];
      }

      array[0] = value;
      array2[0] = value2;
    }
    else
    {
      for ( ; comp(value, array[j - 1]); --j)
      {
        array[j] = array[j - 1];
        array2[j] = array2[j - 1];
      }

      array[j] = value;
      array2[j] = value2;
    }
  }
}

template <typename T, typename T2>
void insertion_sort(T* array, long size, T2* array2)
{
  insertion_sort(array, size, array2, std::less<T>());
}

template <typename T>
void counting_sort(T* array, long size, T maxval)
{
  T nbin = maxval + 1;
  T* count  = (T*) calloc(nbin, sizeof(T));
  T* result = (T*) calloc(size, sizeof(T));
  T* start = (T*) malloc(nbin * sizeof(T));

  // Build histogram of array in count.
  for (long i = 0; i < size; ++i) ++count[array[i]];

  // Put starting location for each bucket in start.
  start[0] = 0;
  for (T i = 1; i < nbin; ++i) start[i] = start[i - 1] + count[i - 1];

  // Put the elements of array into their sorted order in result.
  #pragma mta assert nodep
  for (long i = 0; i < size; ++i)
  {
    T loc = mt_incr(start[array[i]], 1);
    result[loc] = array[i];
  }

  // Copy the sorted array from result to array.
  for (long i = 0; i < size; ++i) array[i] = result[i];

  free(count);
  free(result);
  free(start);
}

template <typename T>
void counting_sort(T* array, unsigned long asize, T maxval)
{
  T nbin = maxval + 1;
  T* count  = (T*) calloc(nbin, sizeof(T));
  T* result = (T*) calloc(asize, sizeof(T));
  T* start = (T*) malloc(nbin * sizeof(T));

  // Build histogram of array in count.
  for (unsigned long i = 0; i < asize; ++i) ++count[array[i]];

  // Put starting location for each bucket in start.
  start[0] = 0;
  for (T i = 1; i < nbin; ++i) start[i] = start[i - 1] + count[i - 1];

  // Put the elements of array into their sorted order in result.
  #pragma mta assert nodep
  for (unsigned long i = 0; i < asize; ++i)
  {
    T loc = mt_incr(start[array[i]], 1);
    result[loc] = array[i];
  }

  // Copy the sorted array from result to array.
  for (unsigned long i = 0; i < asize; ++i) array[i] = result[i];

  free(count);
  free(result);
  free(start);
}

template <typename T, typename T2>
void bucket_sort(T* array, unsigned asize, T maxval, T2* data)
{
  typedef struct ss {
    T key;
    T2 data;
    struct ss* next;
  } bnode;

  int num_b = maxval + 1;
  bnode** buckets = (bnode**) calloc(num_b, sizeof(bnode*));
  T* count = (T*) calloc(num_b, sizeof(T));
  long* b_locks = (long*) calloc(num_b, sizeof(long));
  bnode* nodes = (bnode*) malloc(sizeof(bnode) * asize);
  int pos = 0;

  #pragma mta assert parallel
  for (unsigned i = 0; i < asize; i++)
  {
    int p = mt_incr(pos, 1);
    T key = array[i];
    nodes[p].key = key;
    nodes[p].data = data[i];
    mt_readfe(b_locks[key]);
    nodes[p].next = buckets[key];
    buckets[key] = &nodes[p];
    mt_write(b_locks[key], 1);
    mt_incr(count[key], 1);
  }

  T* start = (T*) malloc(num_b * sizeof(T));
  start[0] = (T) 0;

  for (unsigned i = 1; i < num_b; i++) start[i] = start[i - 1] + count[i - 1];

  T incr = (T) 1;

  #pragma mta assert parallel
  for (unsigned i = 0; i < num_b; i++)
  {
    bnode* tmp = buckets[i];
    int loc = start[i];

    while (tmp)
    {
      array[loc] = tmp->key;
      data[loc] = tmp->data;
      tmp = tmp->next;
      loc = loc + 1;
    }
  }

  free(count);
  free(start);
  free(buckets);
  free(nodes);
  free(b_locks);
}

template <typename T, typename T2>
void bucket_sort_par_cutoff(T* array, unsigned asize, T maxval, T2* data,
                            int par_cutoff)
{
  typedef struct ss {
    T key;
    T2 data;
    struct ss* next;
  } bnode;

  int num_b = maxval + 1;
  bnode** buckets = (bnode**) calloc(num_b, sizeof(bnode*));
  T* count = (T*) calloc(num_b, sizeof(T));
  long* b_locks = (long*) calloc(num_b, sizeof(long));
  bnode* nodes = (bnode*) malloc(sizeof(bnode) * asize);
  int pos = 0;

  if (asize > par_cutoff)
  {
    #pragma mta assert parallel
    for (unsigned i = 0; i < asize; i++)
    {
      int p = mt_incr(pos, 1);
      T key = array[i];
      nodes[p].key = key;
      nodes[p].data = data[i];
      mt_readfe(b_locks[key]);
      nodes[p].next = buckets[key];
      buckets[key] = &nodes[p];
      mt_write(b_locks[key], 1);
      mt_incr(count[key], 1);
    }
  }
  else
  {
    for (unsigned i = 0; i < asize; i++)
    {
      int p = mt_incr(pos, 1);
      T key = array[i];
      nodes[p].key = key;
      nodes[p].data = data[i];
      mt_readfe(b_locks[key]);
      nodes[p].next = buckets[key];
      buckets[key] = &nodes[p];
      mt_write(b_locks[key], 1);
      mt_incr(count[key], 1);
    }
  }

  T* start = (T*) malloc(num_b * sizeof(T));
  start[0] = (T) 0;

  for (unsigned i = 1; i < num_b; i++) start[i] = start[i - 1] + count[i - 1];

  T incr = (T) 1;

  if (num_b > par_cutoff)
  {
    #pragma mta assert parallel
    for (unsigned i = 0; i < num_b; i++)
    {
      bnode* tmp = buckets[i];
      int loc = start[i];

      while (tmp)
      {
        array[loc] = tmp->key;
        data[loc] = tmp->data;
        tmp = tmp->next;
        loc = loc + 1;
      }
    }
  }
  else
  {
    for (unsigned i = 0; i < num_b; i++)
    {
      bnode* tmp = buckets[i];
      int loc = start[i];

      while (tmp)
      {
        array[loc] = tmp->key;
        data[loc] = tmp->data;
        tmp = tmp->next;
        loc = loc + 1;
      }
    }
  }

  free(count);
  free(start);
  free(buckets);
  free(nodes);
  free(b_locks);
}

template <typename T, typename T2>
void bucket_sort(T* array, unsigned asize, T maxval, T2* data, T*& start)
{
  typedef struct ss {
    T key;
    T2 data;
    struct ss* next;
  } bnode;

  unsigned int num_b = maxval + 1;
  bnode** buckets = (bnode**) calloc(num_b, sizeof(bnode*));
  T* count = (T*) calloc(num_b, sizeof(T));
  long* b_locks = (long*) calloc(num_b, sizeof(long));

  bnode* nodes = (bnode*) malloc(sizeof(bnode) * asize);
  long pos = 0;

  #pragma mta assert parallel
  for (unsigned i = 0; i < asize; i++)
  {
    long p = mt_incr(pos, 1);
    T key = array[i];
    nodes[p].key = key;
    nodes[p].data = data[i];
    mt_readfe(b_locks[key]);
    nodes[p].next = buckets[key];
    buckets[key] = &nodes[p];
    mt_write(b_locks[key], 1);
    mt_incr(count[key], 1);
  }

  T* start_local = (T*) malloc(num_b * sizeof(T));
  start_local[0] = (T) 0;

  for (unsigned i = 1; i < num_b; i++)
  {
    start_local[i] =  start_local[i - 1] + count[i - 1];
  }

  start = start_local;
  start_local = 0;

  #pragma mta assert nodep
  for (unsigned i = 0; i < num_b; i++)
  {
    bnode* tmp = buckets[i];
    int loc = start[i];

    while (tmp)
    {
      array[loc] = tmp->key;
      data[loc] = tmp->data;
      tmp = tmp->next;
      loc = loc + 1;
    }
  }

  free(buckets);
  free(nodes);
  free(count);
  free(b_locks);
}

template <typename T, typename T2, typename T3>
void bucket_sort(T* array, unsigned asize, T maxval, T2* data, T3* data2)
{
  typedef struct ss {
    T key;
    T2 data;
    T3 data2;
    struct ss* next;
  } bnode;

  int num_b = maxval + 1;
  bnode** buckets = (bnode**) calloc(num_b, sizeof(bnode*));
  T* count = (T*) calloc(num_b, sizeof(T));
  int* b_locks = (int*) calloc(num_b, sizeof(int));
  bnode* nodes = (bnode*) malloc(sizeof(bnode) * asize);
  int pos = 0;

  #pragma mta assert parallel
  for (unsigned i = 0; i < asize; i++)
  {
    int p = mt_incr(pos, 1);
    T key = array[i];
    nodes[p].key = key;
    nodes[p].data = data[i];
    nodes[p].data2 = data2[i];
    mt_readfe(b_locks[key]);
    nodes[p].next = buckets[key];
    buckets[key] = &nodes[p];
    mt_write(b_locks[key], 1);
    mt_incr(count[key], 1);
  }

  T* start = (T*) malloc(num_b * sizeof(T));
  start[0] = (T) 0;

  for (unsigned i = 1; i < num_b; i++) start[i] = start[i - 1] + count[i - 1];

  T incr = (T) 1;

  #pragma mta assert parallel
  for (unsigned i = 0; i < num_b; i++)
  {
    bnode* tmp = buckets[i];
    int loc = start[i];

    while (tmp)
    {
      array[loc] = tmp->key;
      data[loc] = tmp->data;
      data2[loc] = tmp->data2;
      tmp = tmp->next;
      loc = loc + 1;
    }
  }

  free(count);
  free(start);
  free(buckets);
  free(nodes);
  free(b_locks);
}

template <typename T>
void bucket_sort(T* array, unsigned asize, T maxval)
{
  typedef struct ss {
    T key;
    struct ss* next;
  } bnode;

  int num_b = maxval + 1;
  bnode** buckets = (bnode**) calloc(num_b, sizeof(bnode*));
  T* count = (T*) calloc(num_b, sizeof(T));
  int* b_locks = (int*) calloc(num_b, sizeof(int));
  bnode* nodes = (bnode*) malloc(sizeof(bnode) * asize);
  int pos = 0;

  #pragma mta assert parallel
  for (unsigned i = 0; i < asize; i++)
  {
    int p = mt_incr(pos, 1);
    T key = array[i];
    nodes[p].key = key;
    mt_readfe(b_locks[key]);
    nodes[p].next = buckets[key];
    buckets[key] = &nodes[p];
    mt_write(b_locks[key], 1);
    mt_incr(count[key], 1);
  }

  T* start = (T*) malloc(num_b * sizeof(T));
  start[0] = (T) 0;

  for (unsigned i = 1; i < num_b; i++) start[i] = start[i - 1] + count[i - 1];

  T incr = (T) 1;

  #pragma mta assert parallel
  for (unsigned i = 0; i < num_b; i++)
  {
    bnode* tmp = buckets[i];
    int loc = start[i];

    while (tmp)
    {
      array[loc] = tmp->key;
      tmp = tmp->next;
      loc = loc + 1;
    }
  }

  free(count);
  free(start);
  free(buckets);
  free(nodes);
  free(b_locks);
}

template <typename T, typename T2>
void bucket_sort(T* array, long size, T maxval, T2* data)
{
  typedef struct ss {
    T key;
    T2 data;
    struct ss* next;
  } bnode;

  long num_b = maxval + 1;
  bnode** buckets = (bnode**) calloc(num_b, sizeof(bnode*));
  T* count = (T*) calloc(num_b, sizeof(T));
  long* b_locks = (long*) calloc(num_b, sizeof(long));
  bnode* nodes = (bnode*) malloc(sizeof(bnode) * size);
  long pos = 0;

  #pragma mta assert parallel
  for (long i = 0; i < size; ++i)
  {
    long p = mt_incr(pos, 1);
    T key = array[i];
    nodes[p].key = key;
    nodes[p].data = data[i];
    mt_readfe(b_locks[key]);
    nodes[p].next = buckets[key];
    buckets[key] = &nodes[p];
    mt_write(b_locks[key], 1);
    mt_incr(count[key], 1);
  }

  T* start = (T*) malloc(num_b * sizeof(T));
  start[0] = 0;

  for (long i = 1; i < num_b; ++i) start[i] = start[i - 1] + count[i - 1];

  T incr = 1;

  #pragma mta assert parallel
  for (long i = 0; i < num_b; ++i)
  {
    bnode* tmp = buckets[i];
    T loc = start[i];

    while (tmp)
    {
      array[loc] = tmp->key;
      data[loc] = tmp->data;
      tmp = tmp->next;
      loc = loc + 1;
    }
  }

  free(count);
  free(start);
  free(buckets);
  free(nodes);
  free(b_locks);
}

template <typename T, typename T2>
void bucket_sort_par_cutoff(T* array, long size, T maxval, T2* data,
                            long par_cutoff)
{
  typedef struct ss {
    T key;
    T2 data;
    struct ss* next;
  } bnode;

  long num_b = maxval + 1;
  bnode** buckets = (bnode**) calloc(num_b, sizeof(bnode*));
  T* count = (T*) calloc(num_b, sizeof(T));
  long* b_locks = (long*) calloc(num_b, sizeof(long));
  bnode* nodes = (bnode*) malloc(sizeof(bnode) * size);
  long pos = 0;

  if (size > par_cutoff)
  {
    #pragma mta assert parallel
    for (long i = 0; i < size; ++i)
    {
      long p = mt_incr(pos, 1);
      T key = array[i];
      nodes[p].key = key;
      nodes[p].data = data[i];
      mt_readfe(b_locks[key]);
      nodes[p].next = buckets[key];
      buckets[key] = &nodes[p];
      mt_write(b_locks[key], 1);
      mt_incr(count[key], 1);
    }
  }
  else
  {
    for (long i = 0; i < size; ++i)
    {
      long p = mt_incr(pos, 1);
      T key = array[i];
      nodes[p].key = key;
      nodes[p].data = data[i];
      mt_readfe(b_locks[key]);
      nodes[p].next = buckets[key];
      buckets[key] = &nodes[p];
      mt_write(b_locks[key], 1);
      mt_incr(count[key], 1);
    }
  }

  T* start = (T*) malloc(num_b * sizeof(T));
  start[0] = 0;

  for (long i = 1; i < num_b; ++i) start[i] = start[i - 1] + count[i - 1];

  T incr = 1;

  if (num_b > par_cutoff)
  {
    #pragma mta assert parallel
    for (long i = 0; i < num_b; ++i)
    {
      bnode* tmp = buckets[i];
      T loc = start[i];

      while (tmp)
      {
        array[loc] = tmp->key;
        data[loc] = tmp->data;
        tmp = tmp->next;
        loc = loc + 1;
      }
    }
  }
  else
  {
    for (long i = 0; i < num_b; ++i)
    {
      bnode* tmp = buckets[i];
      T loc = start[i];

      while (tmp)
      {
        array[loc] = tmp->key;
        data[loc] = tmp->data;
        tmp = tmp->next;
        loc = loc + 1;
      }
    }
  }

  free(count);
  free(start);
  free(buckets);
  free(nodes);
  free(b_locks);
}

template <typename T, typename T2, typename T3>
void bucket_sort(T* array, long size, T maxval, T2* data, T3* data2)
{
  typedef struct ss {
    T key;
    T2 data;
    T3 data2;
    struct ss* next;
  } bnode;

  long num_b = maxval + 1;
  bnode** buckets = (bnode**) calloc(num_b, sizeof(bnode*));
  T* count = (T*) calloc(num_b, sizeof(T));
  long* b_locks = (long*) calloc(num_b, sizeof(long));
  bnode* nodes = (bnode*) malloc(sizeof(bnode) * size);
  long pos = 0;

  #pragma mta assert parallel
  for (long i = 0; i < size; ++i)
  {
    long p = mt_incr(pos, 1);
    T key = array[i];
    nodes[p].key = key;
    nodes[p].data = data[i];
    nodes[p].data2 = data2[i];
    mt_readfe(b_locks[key]);
    nodes[p].next = buckets[key];
    buckets[key] = &nodes[p];
    mt_write(b_locks[key], 1);
    mt_incr(count[key], 1);
  }

  T* start = (T*) malloc(num_b * sizeof(T));
  start[0] = 0;

  for (long i = 1; i < num_b; ++i) start[i] = start[i - 1] + count[i - 1];

  T incr = 1;

  #pragma mta assert parallel
  for (long i = 0; i < num_b; ++i)
  {
    bnode* tmp = buckets[i];
    T loc = start[i];

    while (tmp)
    {
      array[loc] = tmp->key;
      data[loc] = tmp->data;
      data2[loc] = tmp->data2;
      tmp = tmp->next;
      loc = loc + 1;
    }
  }

  free(count);
  free(start);
  free(buckets);
  free(nodes);
  free(b_locks);
}

template <typename T>
bool binary_search(T* array, long size, const T& value)
{
  long low = 0;
  long high = size - 1;
  long mid;

  while (low <= high)
  {
    mid = (low + high) / 2;
    if (array[mid] > value)
    {
      high = mid - 1;
    }
    else if (array[mid] < value)
    {
      low = mid + 1;
    }
    else
    {
      return true;
    }
  }

  return false;
}

}

#endif
