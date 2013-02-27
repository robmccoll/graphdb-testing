/*
 *  _________________________________________________________________________
 *  MTGL: The MultiThreaded Graph Library
 *  Copyright (c) 2008 Sandia Corporation.
 *  This software is distributed under the BSD License.
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *  For more information, see the README file in the top MTGL directory.
 *  _________________________________________________________________________
 */

/****************************************************************************/
/*! \file util.hpp

    \brief Common definitions for the mtgl namespace.

    \author Jon Berry (jberry@sandia.gov)
    \author Greg Mackey (gemacke@sandia.gov)
    \author Eric Goodman (elgoodm@sandia.gov)

    \date 12/4/2007
*/
/****************************************************************************/

#ifndef MTGL_UTIL_HPP
#define MTGL_UTIL_HPP

#include <cstdlib>
#include <iostream>

#ifdef __MTA__
#include <sys/mta_task.h>
#include <machine/runtime.h>
#elif defined(_WIN32)
#include <ctime>
#else
#include <ctime>
#include <sys/time.h>
#include <cstring>
#endif

#if defined(USING_QTHREADS)
#include <qthread/qthread.h>
#endif

#include <mtgl/graph_traits.hpp>
#include <mtgl/mtgl_config.h>

// Disable the warning about unknown pragmas.
#ifndef __MTA__
#pragma warning ( disable : 4068 )
#endif

namespace mtgl {

// For specifying algorithm return types.
enum { VERTICES = 0, EDGES };

// Filter types.
const int NO_FILTER =   0;
const int AND_FILTER =  1;
const int PURE_FILTER = 2;

// How edges should be treated.
const int UNDIRECTED = 0;
const int DIRECTED =   1;
const int REVERSED =   2;

/// \brief Reads the value in target, waiting until the "empty bit" is full and
///        setting the "empty bit" to empty.
template <typename T>
T mt_readfe(T& target)
{
#ifdef __MTA__
  return readfe(&target);
#elif USING_QTHREADS
  T ret;
  qthread_readFE(&ret, &target);
  return ret;
#else
  return target;
#endif
}

/// \brief Reads the value in target, waiting until the "empty bit" is full and
///        leaving the "empty bit" full.
template <typename T>
T mt_readff(T& target)
{
#ifdef __MTA__
  return readff(&target);
#elif USING_QTHREADS
  T ret;
  qthread_readFF(&ret, &target);
  return ret;
#else
  return target;
#endif
}

/// \brief Reads the value in target, ignoring the "empty bit".
template <typename T>
T mt_readxx(T& target)
{
#ifdef __MTA__
  return readxx(&target);
#else
  return target;
#endif
}

/// \brief Writes val to target, waiting until the "empty bit" is empty and
///        setting the "empty bit" to full.
template <typename T, typename T2>
void mt_write(T& target, T2 val)
{
#ifdef __MTA__
  writeef(&target, val);
#elif USING_QTHREADS
  qthread_writeEF_const((aligned_t*) &target, (aligned_t) val);
#else
  target = val;
#endif
}


/// \brief Writes val to target, ignoring the "empty bit".
template <typename T, typename T2>
void mt_writexf(T& target, T2 val)
{
#ifdef __MTA__
  writexf(&target, val);
#elif USING_QTHREADS
  qthread_writeF_const((aligned_t*) &target, (aligned_t) val);
#else
  target = val;
#endif
}

/// \brief Performs an atomic integer fetch-and-add.
template <typename T, typename T2>
T mt_incr(T& target, T2 inc)
{
#ifdef __MTA__
  T res = int_fetch_add((int*) &target, inc);
#elif _OPENMP
#ifdef MTGL_ATOMIC_INCR
  T res = __sync_fetch_and_add(&target, inc);
#else
  T res;
  #pragma omp critical
  {
    res = target;
    target += inc;
  }
#endif
#elif USING_QTHREADS
  T res = qthread_incr(&target, inc);
#else
  T res = target;
  target += inc;
#endif

  return res;
}

/// \brief Performs an atomic double fetch-and-add.
template <typename T>
inline double mt_incr(double& target, T inc)
{
#ifdef __MTA__
  double res = mt_readfe(target);
  mt_write(target, res + inc);
#elif _OPENMP
  double res;
  #pragma omp critical
  {
    res = target;
    target += inc;
  }
#elif USING_QTHREADS
  double res = qthread_dincr(&target, inc);
#else
  double res = target;
  target += inc;
#endif

  return res;
}

/// \brief Writes a value of 0 to target and sets its "empty bit" to full.
template <typename T>
void mt_purge(T& target)
{
#ifdef __MTA__
  purge(&target);
#elif USING_QTHREADS
  qthread_empty(&target);
  target = 0;
#else
  target = 0;
#endif
}

class mt_timer {
public:
#ifdef __MTA__
  mt_timer() : ticks(0), freq(mta_clock_freq()) {}
#else
  mt_timer() {}
#endif

#ifdef __MTA__
  void start() { ticks = mta_get_clock(0); }
  void stop() { ticks = mta_get_clock(ticks); }
  double getElapsedSeconds() { return ticks / freq; }
  long getElapsedTicks() { return ticks; }

#elif defined(_WIN32)
  void start() { start_time = std::clock(); }
  void stop() { stop_time = std::clock(); }

  double getElapsedSeconds()
  {
    std::clock_t ticks = stop_time - start_time;
    return ticks / (double) CLOCKS_PER_SEC;
  }

  std::clock_t getElapsedTicks()
  {
    std::clock_t ticks = stop_time - start_time;
    return ticks;
  }

#else
  void start() { gettimeofday(&start_time, NULL); }
  void stop() { gettimeofday(&stop_time, NULL); }

  long getElapsedTicks()
  {
    return (long) (getElapsedSeconds() * CLOCKS_PER_SEC);
  }

  double getElapsedSeconds()
  {
    double start_t = start_time.tv_sec + start_time.tv_usec * 1e-6;
    double stop_t = stop_time.tv_sec + stop_time.tv_usec * 1e-6;
    return stop_t - start_t;
  }
#endif

private:
#ifdef __MTA__
  int start_time;
  int stop_time;
  int ticks;
  double freq;
#elif defined(_WIN32)
  std::clock_t start_time;
  std::clock_t stop_time;
#else
  struct timeval start_time;
  struct timeval stop_time;
#endif
};

template <typename T>
void init_mta_counters(mt_timer& timer, T& issues, T& memrefs, T& concur,
                       T& streams)
{
#ifdef __MTA__
  issues = mta_get_task_counter(RT_ISSUES);
  memrefs = mta_get_task_counter(RT_MEMREFS);
  concur = mta_get_task_counter(RT_CONCURRENCY);
  streams = mta_get_task_counter(RT_STREAMS);
#endif
  timer.start();
}

template <typename T>
void init_mta_counters(mt_timer& timer, T& issues, T& memrefs, T& concur,
                       T& streams, T& traps)
{
#ifdef __MTA__
  issues = mta_get_task_counter(RT_ISSUES);
  memrefs = mta_get_task_counter(RT_MEMREFS);
  concur = mta_get_task_counter(RT_CONCURRENCY);
  streams = mta_get_task_counter(RT_STREAMS);
  traps = mta_get_task_counter(RT_TRAP);
#endif
  timer.start();
}

template <typename T>
void sample_mta_counters(mt_timer& timer, T& issues, T& memrefs, T& concur,
                         T& streams)
{
#ifdef __MTA__
  issues = mta_get_task_counter(RT_ISSUES) - issues;
  memrefs = mta_get_task_counter(RT_MEMREFS) - memrefs;
  concur = mta_get_task_counter(RT_CONCURRENCY) - concur;
  streams = mta_get_task_counter(RT_STREAMS) - streams;
#endif
  timer.stop();
}

template <typename T>
void sample_mta_counters(mt_timer& timer, T& issues, T& memrefs, T& concur,
                         T& streams, T& traps)
{
#ifdef __MTA__
  issues = mta_get_task_counter(RT_ISSUES) - issues;
  memrefs = mta_get_task_counter(RT_MEMREFS) - memrefs;
  concur = mta_get_task_counter(RT_CONCURRENCY) - concur;
  streams = mta_get_task_counter(RT_STREAMS) - streams;
  traps = mta_get_task_counter(RT_TRAP) - traps;
#endif
  timer.stop();
}

template <typename T>
void print_mta_counters(mt_timer& timer, int m,
                        T& issues, T& memrefs, T& concur, T& streams)
{
#ifdef __MTA__
  std::cout << "secs: " << timer.getElapsedSeconds() << ", issues: "
            << issues << ", memrefs: " << memrefs << ", concurrency: "
            << concur / static_cast<double>(timer.getElapsedTicks())
            << ", streams: "
            << streams / static_cast<double>(timer.getElapsedTicks())
            << std::endl
            << "memrefs/edge: " << memrefs / static_cast<double>(m)
            << std::endl;
#else
  std::cout << "secs: " << timer.getElapsedSeconds() << std::endl;
#endif
}

template <typename T>
void print_mta_counters(mt_timer& timer, int m,
                        T& issues, T& memrefs, T& concur, T& streams, T& traps)
{
#ifdef __MTA__
  std::cout << "secs: " << timer.getElapsedSeconds() << ", issues: "
            << issues << ", memrefs: " << memrefs << ", concurrency: "
            << concur / static_cast<double>(timer.getElapsedTicks())
            << ", streams: "
            << streams / static_cast<double>(timer.getElapsedTicks())
            << ", traps: " << traps << std::endl
            << "memrefs/edge: " << memrefs / static_cast<double>(m)
            << std::endl;
#else
  std::cout << "secs: " << timer.getElapsedSeconds() << std::endl;
#endif
}

// Note: We can't include <utility> to get std::pair because of MTA issues;
//       instead, we duplicate in our own namespace.
template <typename T1, typename T2>
class pair {
public:
  typedef T1 first_type;
  typedef T2 second_type;

  pair() : first(T1()), second(T2()) {}
  pair(const T1& x, const T2& y) : first(x), second(y) {}
  pair(const pair<T1, T2>& p) : first(p.first), second(p.second) {}

  template <typename U, typename V>
  pair(const pair<U, V>& p) : first(p.first), second(p.second) {}

  int operator==(const pair<T1, T2>& p)  const
  {
    return (first == p.first && second == p.second);
  }

  friend std::ostream& operator<<(std::ostream& os, const pair& p)
  {
    os << "[" << p.first << "," << p.second << "]";

    return os;
  }

public:
  T1 first;
  T2 second;
};

template <typename T1, typename T2>
class pair<T1&, T2&> {
public:
  typedef T1 first_type;
  typedef T2 second_type;

  pair() : first(T1()), second(T2()) {}
  pair(T1& x, T2& y) : first(x), second(y) {}
  pair(pair<T1, T2>& p) : first(p.first), second(p.second) {}

  template <typename U, typename V>
  pair(pair<U, V>& p) : first(p.first), second(p.second) {}

  pair<T1&, T2&>& operator=(const pair<T1&, T2&>& a)
  {
    first = a.first;
    second = a.second;
    return *this;
  }

  pair<T1&, T2&>& operator=(const pair<T1, T2>& a)
  {
    first = a.first;
    second = a.second;
    return *this;
  }

  int operator==(pair<T1, T2>& p)  const
  {
    return (first == p.first && second == p.second);
  }

  friend std::ostream& operator<<(std::ostream& os, const pair& p)
  {
    os << "[" << p.first << "," << p.second << "]";

    return os;
  }

public:
  T1& first;
  T2& second;
};

template <typename T1, typename T2, typename T3>
class triple {
public:
  typedef T1 first_type;
  typedef T2 second_type;
  typedef T3 third_type;

  triple() : first(T1()), second(T2()), third(T3()) {}

  triple(const T1& x, const T2& y, const T3& z) :
    first(x), second(y), third(z) {}

  triple(const triple<T1, T2, T3>& p) :
    first(p.first), second(p.second), third(p.third) {}

  template <typename U, typename V, typename W>
  triple(const triple<U, V, W>& p) :
    first(p.first), second(p.second), third(p.third) {}

  int operator==(const triple<T1, T2, T3>& p)  const
  {
    return (first == p.first && second == p.second && third == p.third);
  }

  bool operator<(const triple<T1, T2, T3>& p) const
  {
    return first < p.first ||
           (!(p.first < first) && (second < p.second ||
                                   (!(p.second < second) && third < p.third)));
  }

  friend std::ostream& operator<<(std::ostream& os, const triple& p)
  {
    os << "[" << p.first << "," << p.second << "," << p.third << "]";

    return os;
  }

public:
  T1 first;
  T2 second;
  T3 third;
};

template <typename T1, typename T2, typename T3>
class triple<T1&, T2&, T3&> {
public:
  typedef T1 first_type;
  typedef T2 second_type;
  typedef T3 third_type;

  triple() : first(T1()), second(T2()), third(T3()) {}

  triple(T1& x, T2& y, T3& z) : first(x), second(y), third(z) {}

  triple(triple<T1, T2, T3>& p) :
    first(p.first), second(p.second), third(p.third) {}

  template <typename U, typename V, typename W>
  triple(triple<U, V, W>& p) :
    first(p.first), second(p.second), third(p.third) {}

  triple<T1&, T2&, T3&>& operator=(const triple<T1&, T2&, T3&>& a)
  {
    first = a.first;
    second = a.second;
    third = a.third;
    return *this;
  }

  triple<T1&, T2&, T3&>& operator=(const triple<T1, T2, T3>& a)
  {
    first = a.first;
    second = a.second;
    third = a.third;
    return *this;
  }

  int operator==(triple<T1, T2, T3>& p)  const
  {
    return (first == p.first && second == p.second && third == p.third);
  }

  friend std::ostream& operator<<(std::ostream& os, const triple& p)
  {
    os << "[" << p.first << "," << p.second << "," << p.third << "]";

    return os;
  }

public:
  T1& first;
  T2& second;
  T3& third;
};

template <typename T1, typename T2, typename T3, typename T4>
class quadruple {
public:
  typedef T1 first_type;
  typedef T2 second_type;
  typedef T3 third_type;
  typedef T4 fourth_type;

  quadruple() : first(T1()), second(T2()), third(T3()), fourth(T4()) {}

  quadruple(const T1& w, const T2& x, const T3& y, const T4& z) :
    first(w), second(x), third(y), fourth(z) {}

  quadruple(const quadruple<T1, T2, T3, T4>& q) :
    first(q.first), second(q.second), third(q.third), fourth(q.fourth) {}

  template <typename U, typename V, typename W, typename X>
  quadruple(const quadruple<U, V, W, X>& q) :
    first(q.first), second(q.second), third(q.third), fourth(q.fourth) {}

  int operator==(const quadruple<T1, T2, T3, T4>& q) const
  {
    return (first == q.first && second == q.second && third == q.third &&
            fourth == q.fourth);
  }

  friend std::ostream& operator<<(std::ostream& os, const quadruple& q)
  {
    os << "[" << q.first << "," << q.second << "," << q.third << ","
        << q.fourth << "]";

    return os;
  }

public:
  T1 first;
  T2 second;
  T3 third;
  T4 fourth;
};

template <typename T1, typename T2, typename T3, typename T4>
class quadruple<T1&, T2&, T3&, T4&> {
public:
  typedef T1 first_type;
  typedef T2 second_type;
  typedef T3 third_type;
  typedef T4 fourth_type;

  quadruple() : first(T1()), second(T2()), third(T3()), fourth(T4()) {}

  quadruple(T1& w, T2& x, T3& y, T4& z) :
    first(w), second(x), third(y), fourth(z) {}

  quadruple(quadruple<T1, T2, T3, T4>& p) :
    first(p.first), second(p.second), third(p.third), fourth(p.fourth) {}

  template <typename U, typename V, typename W, typename X>
  quadruple(quadruple<U, V, W, X>& p) :
    first(p.first), second(p.second), third(p.third), fourth(p.fourth) {}

  quadruple<T1&, T2&, T3&, T4&>&
  operator=(const quadruple<T1&, T2&, T3&, T4&>& a)
  {
    first = a.first;
    second = a.second;
    third = a.third;
    fourth = a.fourth;
    return *this;
  }

  quadruple<T1&, T2&, T3&, T4&>& operator=(const quadruple<T1, T2, T3, T4>& a)
  {
    first = a.first;
    second = a.second;
    third = a.third;
    fourth = a.fourth;
    return *this;
  }

  int operator==(quadruple<T1, T2, T3, T4>& p)  const
  {
    return (first == p.first && second == p.second && third == p.third &&
            fourth == p.fourth);
  }

  friend std::ostream& operator<<(std::ostream& os, const quadruple& p)
  {
    os << "[" << p.first << "," << p.second << "," << p.third << ","
        << p.fourth << "]";

    return os;
  }

public:
  T1& first;
  T2& second;
  T3& third;
  T4& fourth;
};

template <typename T1, typename T2>
inline pair<T1&, T2&> tie(T1& t1, T2& t2)
{
  return pair<T1&, T2&> (t1, t2);
}

template <typename T1, typename T2, typename T3>
inline triple<T1&, T2&, T3&> tie(T1& t1, T2& t2, T3& t3)
{
  return triple<T1&, T2&, T3&> (t1, t2, t3);
}

template <typename T1, typename T2, typename T3, typename T4>
inline quadruple<T1&, T2&, T3&, T4&> tie(T1& t1, T2& t2, T3& t3, T4& t4)
{
  return quadruple<T1&, T2&, T3&, T4&> (t1, t2, t3, t4);
}

template <typename T>
static void order_pair(T& a, T& b)
{
  if (a > b)
  {
    T tmp = a;
    a = b;
    b = tmp;
  }
}

/*! \fn record_time(mt_timer& timer, double& total_time, const char* message)

    \brief Found myself (ELG) writing this same bit of code a lot, so
           put it into a function.  It stops a timer, gets the time,
           and then outputs a message.

    \param timer The timer used to get the time since the last start().
    \param message The message to output to stdout.
    \return Returns the elapsed seconds.

    \author Eric Goodman
    \date 2/2011
*/
inline double record_time(mt_timer& timer, const char* message)
{
  #pragma mta fence
  timer.stop();

  double t = timer.getElapsedSeconds();
  std::cout << message << t << std::endl;

  return t;
}

}

#include <mtgl/snap_util.h>

#endif
