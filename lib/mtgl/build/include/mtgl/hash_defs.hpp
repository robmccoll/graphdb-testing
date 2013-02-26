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
/*! \file hash_defs.hpp

    \brief Definition of default hash functions, key equality functions,
           and update function objects for hash tables.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 3/31/2010
*/
/****************************************************************************/

#ifndef MTGL_HASH_DEFS_HPP
#define MTGL_HASH_DEFS_HPP

#include <string>

#include <mtgl/types.hpp>
#include <mtgl/util.hpp>

namespace mtgl {

/// \brief  Function object that serves as the default equality comparison
///         for hash sets and hash tables.
template <typename T>
struct default_eqfcn {
  bool operator()(const T& t1, const T& t2) const { return t1 == t2; }
};

template <>
struct default_eqfcn<char*> {
  bool operator()(const char* t1, const char* t2) const
  { return strcmp(t1, t2) == 0; }
};

template <>
struct default_eqfcn<const char*> {
  bool operator()(const char* t1, const char* t2) const
  { return strcmp(t1, t2) == 0; }
};

/// Definition for size_type used by hash sets and hash tables.
typedef unsigned long hash_size_type;

/// Default hash function for strings, char arrays, etc.
template <typename T>
inline
hash_size_type string_hash_func(const char* key, const T& key_len)
{
  hash_size_type hash = 0;

  for (T i = 0; i < key_len; ++i)
  {
    hash = key[i] + (hash << 6) + (hash << 16) - hash;
  }

  return hash;
}

/// Default hash function for all integer types.
template <typename T>
inline
hash_size_type integer_hash_func(const T& key)
{
/*
    hash_size_type k = static_cast<hash_size_type>(key);

    k = (~k) + (k << 21); // k = (k << 21) - k - 1;
    k = k ^ (k >> 24);
    k = (k + (k << 3)) + (k << 8); // k * 265
    k = k ^ (k >> 14);
    k = (k + (k << 2)) + (k << 4); // k * 21
    k = k ^ (k >> 28);
    k = k + (k << 31);

    return k;
*/

//    return static_cast<hash_size_type>(key);
#if SIZEOF_UNSIGNED_LONG == 8
  return static_cast<hash_size_type>(key) * 31280644937747LL;
#else
  return static_cast<hash_size_type>(key) * 19922923;
#endif
}

/// \brief  Function object that serves as the default hash function for a
///         variety of key types for hash_sets and hash_tables.
template <typename T>
struct default_hash_func {};

template <>
struct default_hash_func<short>
{
  hash_size_type operator()(short key) const
  { return integer_hash_func(key); }
};

template <>
struct default_hash_func<unsigned short>
{
  hash_size_type operator()(unsigned short key) const
  { return integer_hash_func(key); }
};

template <>
struct default_hash_func<int>
{
  hash_size_type operator()(int key) const
  { return integer_hash_func(key); }
};

template <>
struct default_hash_func<unsigned int>
{
  hash_size_type operator()(unsigned int key) const
  { return integer_hash_func(key); }
};

template <>
struct default_hash_func<long>
{
  hash_size_type operator()(long key) const
  { return integer_hash_func(key); }
};

template <>
struct default_hash_func<unsigned long>
{
  hash_size_type operator()(unsigned long key) const
  { return integer_hash_func(key); }
};

template <>
struct default_hash_func<long long>
{
  hash_size_type operator()(long long key) const
  { return integer_hash_func(key); }
};

template <>
struct default_hash_func<unsigned long long>
{
  hash_size_type operator()(unsigned long long key) const
  { return integer_hash_func(key); }
};

template <>
struct default_hash_func<std::string> {
public:
  hash_size_type operator()(const std::string& key) const
  { return string_hash_func(key.c_str(), key.size()); }
};

template <>
struct default_hash_func<const std::string> {
public:
  hash_size_type operator()(const std::string& key) const
  { return string_hash_func(key.c_str(), key.size()); }
};

template <>
struct default_hash_func<char*> {
public:
  hash_size_type operator()(const char* key) const
  { return string_hash_func(key, strlen(key)); }
};

template <>
struct default_hash_func<const char*> {
public:
  hash_size_type operator()(const char* key) const
  { return string_hash_func(key, strlen(key)); }
};

template <>
struct default_hash_func<char> {
public:
  hash_size_type operator()(char key) const
  { return integer_hash_func(key); }
};

template <>
struct default_hash_func<unsigned char> {
public:
  hash_size_type operator()(unsigned char key) const
  { return integer_hash_func(key); }
};

template <>
struct default_hash_func<signed char> {
public:
  hash_size_type operator()(signed char key) const
  { return integer_hash_func(key); }
};

/// Function object for update functions that performs an mt_incr.
template <typename T>
struct hash_mt_incr
{
  void operator()(T& v, const T& value) const
  { mt_incr(v, value); }
};

/// \brief Function object for update functions that updates a value only if
///        the new value is less than the current value.
template <typename T>
struct hash_min
{
  void operator()(T& v, const T& value) const
  {
    if (value < v)
    {
      T probe = mt_readfe(v);
      if (value < v)
      {
        mt_write(v, value);
      }
      else
      {
        mt_write(v, probe);
      }
    }
  }
};

/// \brief Assigns contiguous ids to the values for the entries in the table.
///        returns the number of keys.
///
/// \param ht The hash table.
/// \param starting_value The starting value for the contiguous numbering.
template <typename HT>
inline
void assign_contiguous_ids(HT& ht, typename HT::size_type starting_value = 0)
{
  typename HT::size_type stream_id = 0;
  typename HT::size_type num_streams = 1;

  #pragma mta use 75 streams
  #pragma mta for all streams stream_id of num_streams
  {
    typename HT::thread_iterator begin_iter = ht.begin(stream_id, num_streams);
    typename HT::thread_iterator end_iter = ht.end(stream_id, num_streams);

    // Claim a block of indices for the occupied entries in the stream's
    // region.
    typename HT::size_type start_index = mt_incr(starting_value,
                                                 end_iter - begin_iter);

    // Assign contiguous ids to the occupied entries in this stream's
    // region.
    for ( ; begin_iter != end_iter; ++begin_iter)
    {
      begin_iter->second = start_index;
      ++start_index;
    }
  }
}

/// \brief Puts the keys from the hash table in a contiguous array, and
///        returns the number of keys.
///
/// \param ht The hash table.
/// \param keys The keys array that gets overwritten.  This is assumed to be
///             size() or bigger.
/// \return The number of keys.
template <typename HT>
inline
typename HT::size_type
get_keys(HT& ht, typename HT::key_type* keys)
{
  typename HT::size_type num_keys = 0;

  typename HT::size_type stream_id = 0;
  typename HT::size_type num_streams = 1;

  #pragma mta use 100 streams
  #pragma mta for all streams stream_id of num_streams
  {
    typename HT::thread_iterator begin = ht.begin(stream_id, num_streams);
    typename HT::thread_iterator end = ht.end(stream_id, num_streams);

    // Claim a block of indices for the occupied entries in the stream's
    // region.
    typename HT::size_type start_index = mt_incr(num_keys, end - begin);

    // Get the keys from this stream's region.
    for ( ; begin != end; ++begin)
    {
      keys[start_index] = begin->first;
      ++start_index;
    }
  }

  return num_keys;
}

// Algorithm: P.J. Weinberg
// Impl: M. Neumann, 1989
// Modified: J. Berry, 2011
template <typename size_type>
size_type hashpjw(size_type M, char* t, int numbytes)
{
  uint64_t h = 0;
  uint64_t g;

  for (int i = 0; i < numbytes; ++i, ++t)
  {
    h = (h << 4) + *t;
    g = h & 0xf0000000;

    if (g)
    {
      h ^= g >> 24;
      h ^= g;
    }
  }

  return ((size_type) h ) % M;
}

// Algorithm: FowlerNollVo
// Impl: J. Berry, 2011 (Based on wikipedia entry)
template <typename size_type>
size_type fnv1(size_type t)
{
  printf("sizeof(size_type): %lu\n", sizeof(size_type));
  fflush(stdout);
  assert(sizeof(size_type) == 8);

  size_type hash = 14695981039346656037UL;
  uint8_t* octets = (uint8_t*) &t;

  for (int i = 0; i < 8; ++i)
  {
    uint8_t octet = octets[i];
    hash = hash * 1099511628211UL;
    hash = hash ^ octet;
  }

  return hash;
}

}

#endif
