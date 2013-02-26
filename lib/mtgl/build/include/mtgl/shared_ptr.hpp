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
/*! \file shared_ptr.hpp

    \brief Implementation of a reference counting pointer.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 12/7/2010
*/
/****************************************************************************/

#ifndef MTGL_SHARED_PTR_HPP
#define MTGL_SHARED_PTR_HPP

#include <mtgl/util.hpp>

namespace mtgl {

namespace detail {

struct dynamic_cast_tag {};

}

template <typename T>
class shared_ptr {
public:
  typedef T element_type;

  template <typename Y> friend class shared_ptr;

  explicit shared_ptr(T* p = 0) : ptr(p), counter(0)
  {
    if (ptr != 0)
    {
      counter = new unsigned long(1);
    }
  }

  shared_ptr(const shared_ptr& p) : ptr(p.ptr), counter(p.counter)
  {
    if (ptr != 0) mt_incr(*counter, 1);
  }

  template <typename Y>
  shared_ptr(const shared_ptr<Y>& p) : ptr(static_cast<T*>(p.ptr)),
                                       counter(p.counter)
  {
    if (ptr != 0) mt_incr(*counter, 1);
  }

  template <typename Y>
  shared_ptr(shared_ptr<Y> const & p, detail::dynamic_cast_tag):
    ptr(dynamic_cast<T*>(p.ptr)), counter(p.counter)
  {
    if (ptr == 0) // need to allocate new counter -- the cast failed
    {
      counter = new unsigned long(1);
    }
    else
    {
      mt_incr(*counter, 1);
    }
  }

  ~shared_ptr()
  {
    if (ptr != 0 && mt_incr(*counter, -1) == 1)
    {
      delete ptr;
      delete counter;
    }
  }

  shared_ptr& operator=(const shared_ptr& rhs)
  {
    if (ptr != 0 && mt_incr(*counter, -1) == 1)
    {
      delete ptr;
      delete counter;
    }

    ptr = rhs.ptr;
    counter = rhs.counter;

    if (ptr != 0) mt_incr(*counter, 1);

    return *this;
  }

  T& operator*() const { return *ptr; }
  T* operator->() const { return ptr; }
  T* get() const { return ptr; }

  operator bool () const { return ptr != 0; }

  void reset()
  {
    if (ptr != 0 && mt_incr(*counter, -1) == 1)
    {
      delete ptr;
      delete counter;
    }

    ptr = 0;
    counter = 0;
  }

  template <typename Y>
  void reset(Y* p)
  {
    // Don't reset if the new pointer is the same as the current.
    if (ptr == 0 || ptr != p)
    {
      if (ptr != 0 && mt_incr(*counter, -1) == 1)
      {
        delete ptr;
        delete counter;
      }

      ptr = static_cast<T*>(p);
      counter = new unsigned long(1);
    }
  }

  bool unique() const { return ptr != 0 && *counter == 1; }

  friend inline bool operator<(const shared_ptr& a, const shared_ptr& b)
  { return a.get() < b.get(); }

  friend inline bool operator==(const shared_ptr& a, const shared_ptr& b)
  { return a.get() == b.get(); }

  friend inline bool operator!=(const shared_ptr& a, const shared_ptr& b)
  { return a.get() != b.get(); }

private:
  T* ptr;
  unsigned long* counter;
};

template<typename T, typename U>
inline
shared_ptr<T> dynamic_pointer_cast(shared_ptr<U> const &p)
{
  return shared_ptr<T>(p, detail::dynamic_cast_tag());
}

}

#endif
