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
/*! \file mmap_traits.hpp

    \brief This simple class implements a traits class for mmapping data
           that gives the type of data to be mmapped.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 2/14/2011
*/
/****************************************************************************/

#ifndef MTGL_MMAP_TRAITS_HPP
#define MTGL_MMAP_TRAITS_HPP

namespace mtgl {

const unsigned long MMAP_TYPE_NOT_DEFINED = 0;

template <typename T>
class mmap_traits {
public:
  static const unsigned long type = MMAP_TYPE_NOT_DEFINED;
};

template <>
class mmap_traits<bool> {
public:
  static const unsigned long type = 1;
};

template <>
class mmap_traits<char> {
public:
  static const unsigned long type = 2;
};

template <>
class mmap_traits<signed char> {
public:
  static const unsigned long type = 3;
};

template <>
class mmap_traits<unsigned char> {
public:
  static const unsigned long type = 4;
};

template <>
class mmap_traits<short> {
public:
  static const unsigned long type = 5;
};

template <>
class mmap_traits<unsigned short> {
public:
  static const unsigned long type = 6;
};

template <>
class mmap_traits<int> {
public:
  static const unsigned long type = 7;
};

template <>
class mmap_traits<unsigned> {
public:
  static const unsigned long type = 8;
};

template <>
class mmap_traits<long> {
public:
  static const unsigned long type = 9;
};

template <>
class mmap_traits<unsigned long> {
public:
  static const unsigned long type = 10;
};

template <>
class mmap_traits<long long> {
public:
  static const unsigned long type = 11;
};

template <>
class mmap_traits<unsigned long long> {
public:
  static const unsigned long type = 12;
};

template <>
class mmap_traits<float> {
public:
  static const unsigned long type = 13;
};

template <>
class mmap_traits<double> {
public:
  static const unsigned long type = 14;
};

template <>
class mmap_traits<long double> {
public:
  static const unsigned long type = 15;
};

}

#endif
