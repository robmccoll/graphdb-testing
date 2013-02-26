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
/*! \file qalloc.h

    \brief Basic statistics.

    \author Jon Berry (jberry@sandia.gov)

    \date 6/2/2008
*/
/****************************************************************************/

#ifndef MTGL_STATS_HPP
#define MTGL_STATS_HPP

#include <cmath>

namespace mtgl {

template <typename T>
double E(T* a, int sz)
{
  double total = 0.0;

  #pragma mta assert nodep
  for (int i = 0; i < sz; i++) total += (double) a[i];

  return total / sz;
}

/*!
  Harmonic mean = sz / (1/a[0] + 1/a[1] + ... + 1/a[sz-1])
  \author Eric L. Goodman (elgoodm@sandia.gov)
 */
template <typename T>
double harmonic_E(T* a, int sz)
{
  double total = 0.0;

  #pragma mta assert nodep
  for (int i = 0; i < sz; i++) total += (1 / (double) a[i]);

  return sz / total;
}

template <typename T>
double var(T* a, int sz)
{
  double sum_sqr = 0.0;
  double exp_val = E(a, sz);

  #pragma mta assert nodep
  for (int i = 0; i < sz; i++) sum_sqr += ((double) a[i] * (double) a[i]);

  return (sum_sqr / (double) sz) - (exp_val * exp_val);
}

/*!
  Harmonic variant of variance.
  \author Eric L. Goodman (elgoodm@sandia.gov)
 */
template <typename T>
double harmonic_var(T* a, int sz)
{
  double sum_sqr = 0.0;
  double exp_val = harmonic_E(a, sz);

  #pragma mta assert nodep
  for (int i = 0; i < sz; i++) sum_sqr += ((double) a[i] * (double) a[i]);

  return (sum_sqr / (double) sz) - (exp_val * exp_val);
}

template <typename T>
double std_dev(T* a, int sz)
{
  return sqrt(var(a, sz));
}

/*!
  Harmonic variant of standard deviation.
  \author Eric L. Goodman (elgoodm@sandia.gov)
 */
template <typename T>
double harmonic_std_dev(T* a, int sz)
{
  return sqrt(harmonic_var(a, sz));
}



template <typename T>
double correlation_coef(T* a, T* b, int sz)
{
  double e_a = E(a, sz);
  double e_b = E(b, sz);
  double* ab = new double[sz];

  #pragma mta assert nodep
  for (int i = 0; i < sz; i++) ab[i] = (double) a[i] * (double) b[i];

  double e_ab = E(ab, sz);
  double stda = std_dev(a, sz);
  double stdb = std_dev(b, sz);
  double denom = stda * stdb;

  delete [] ab;

  return (e_ab - (e_a * e_b)) / denom;
}

/*!
  Can be used to get quartiles, median, and general percentiles.
  Assumes a sorted array in ascending order.
  \param a The sorted array.
  \param sz The length of the array.
  \param percentile e.g. 0.5 for median, 0.25 for first quartile, etc.
  \author Eric L. Goodman (elgoodm@sandia.gov)
 */
template <typename T, typename size_type>
double percentile(T* a, size_type sz, double percentile) 
{
  double index = percentile * sz;
  if (floor(index) - index != 0) {
    unsigned long i = floor(index);
    if (i > 0) {
      return a[i-1];
    } else {
      return a[0];
    }
  } else {
    unsigned long i = static_cast<long>(index);
    if (i > 0) {
      return static_cast<double>((a[i] + a[i-1])) / 2;
    } else {
      return a[0];
    }
  }
}

}

#endif
