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

// *****************************************************************
// Algorithm Copyright Pacific Northwest National Laboratories 2009.
// *****************************************************************

/****************************************************************************/
/*! \file hachar.hpp

    \brief A map data structure implemented using hashing with chaining and
           region-based memory allocation.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 7/14/2010
*/
/****************************************************************************/

#ifndef HACHAR_HPP
#define HACHAR_HPP

#include <sstream>

#include <mtgl/util.hpp>
#include <mtgl/partitioning.hpp>
#include <mtgl/hash_defs.hpp>

namespace mtgl {

/*! \brief A map data structure implemented using hashing with chaining and
           region-based memory allocation.

    \author Greg Mackey (gemacke@sandia.gov)

    \tparam K Key type.
    \tparam t Value type.
    \tparam HF Function object type defining the hash function.
    \tparam EQF Function object type defining the equality function for keys.

    This hash table is based upon the HACHAR hash table described in the paper
    "Hashing Strategies for the Cray XMT" by Goodman, et. al.  The
    algorithm was primarily developed by David Haglin at Pacific Northwest
    National Labs.  The hash table uses fine-grained locking at the entry
    level, and it uses the double locking method described in the paper.  It
    also avoids locking whenever possible.  Locking is only used when
    modifying a key.  All of these methods increase the parallelism of the
    hash table.

    The hash table enforces unique keys.  Note that keys cannot be modified.
    Once they are inserted into the hash, they are there for the lifetime of
    the table.  The values, however, are modifiable via the update or the
    visit functions.  Because performing an int_fetch_add during insert on an
    element counter causes a hotspot, the table doesn't keep track of the
    number of elements currently in the table.  Instead, a linear cost size()
    function is provided.

    The [] operator is slightly less scalable and slightly slower than the
    other element access functions (insert, update, etc.).  Those functions
    should be preferred over the [] operator when possible.

    There are two main classes of functions: those designed to be called in a
    parallel context and those designed to be called in a serial context.  The
    parallel context functions can be further divided into functions that
    modify the keys and functions that don't.  Here is the division of
    functions:

      - Parallel Context Functions
         - Modify the Key
            - operator[]()
            - insert()
            - update_insert()

         - Don't Modify the Key
            - member()
            - lookup()
            - update()

      - Serial Context Functions
         - Constructors
         - Destructor
         - operator=()
         - size()
         - empty()
         - swap()
         - clear()
         - visit()
         - print()
         - print_stats()
         - print_array_chunk_structure()

    Parallel context functions that don't modify the key are sometimes
    parallelized automatically by the compiler.  Sometimes, though, they need
    a "#pragma mta assert nodep" line before the loop containing the
    function.  Parallel context functions the modify the key require a
    "#pragma mta assert parallel" statement before the loop containing the
    function.  This is because these functions implement locking.  Most
    functions seem to get the best performance when block scheduled.  Here is
    an example of parallelizing insert, a parallel context function that
    modifies the key.

      #pragma mta block schedule
      #pragma mta assert parallel
      for (int i = 0; i < 1000; ++i)
      {
        xht.insert(i, i);
      }

    Here is an example of parallelizing update, a parallel context function
    that doesn't modify the key.

      #pragma mta block schedule
      #pragma mta assert nodep
      for (int i = 0; i < 1000; ++i)
      {
        xht.update(i, i * 2);
      }

    All of the parallel context functions can be interleaved with each other
    inside loops.  They are guaranteed to successfully complete to the
    expected state even though locking is only used around modifying the key
    and not around modifying the value.

    One thing to be careful of is if two threads update the same entry inside
    a single parallel region.  Without some form of locking, you could get
    trash in the value of the entry.  However, this is the normal behavior on
    the XMT whenever two threads update any single memory location in the same
    parallel region.  With the [] operator, you can use the locking primitives
    directly as in:

      mt_incr(ht[3], 2);

    You can provide locking with update() if you use the form, described
    below, that takes a visitor.  The locking would occur in the visitor.

    We provide versions of the update functions that take a third parameter
    of a function object that is applied to the value.  The second parameter
    for the update function is passed to the function object.  This provides
    a lot of flexibility for the user.  There are two pre-defined function
    objects for use with update in hash_defs.hpp: hash_mt_incr and hash_min.
    The functor hash_mt_incr peforms an atomic int-fetch-add on the value
    incrementing it by the second parameter passed to update.  The functor
    hash_min updates value with the second parameter passed to update only if
    the second parameter is less than the current value.  Here is the
    definition of hash_min as an example update function object.  All update
    function objects need to define the operator() as the example below where
    'v' is the current value in the hash table and 'value' is the second
    parameter to update.  Note that whatever locking is necessary on v to
    provide accurate results must be implemented by the function object.

      template <typename T>
      struct hash_mt_incr
      {
        void operator()(T& v, const T& value) const
        { mt_incr(v, value); }
      };

    The hash table provides two types of iteration: loop-based and
    stream-based.  Loop-based iteration is provided by the visit function.
    Visitors are used for loop-based parallelism because iterators for single
    element access cause considerable complication with thread safety.
    The visit function loops over the table and applies a visitor function
    object to each occupied element in the table.  The visitor has access
    to a single element at a time.  The () operator of the visitor should
    have a single parameter that is a pointer to a value_type.  The value_type
    is a pair of the key and the value.  The key is the 'first' field, and the
    value is the second field.  Below is an example of a visitor object that
    increments the value by the key.

      template<typename HT>
      class table_visitor {
      public:
        void operator()(typename HT::value_type* i) { i->second += i->first; }
      };

    Now an example of using the visit function.

      hachar<int, int> ht;
      <Items inserted into the hash here.>
      table_visitor<hachar<int, int> > tv;
      ht.visit(tv);

    Stream-based iteration is providing using thread_iterators in conjunction
    with the "for all streams i of n" pragma.  Since the code inside a "for
    all streams" region is executed serially by each stream, the
    thread_iterators can hold state information. The user gets
    thread_iterators by using the begin and end functions.  The begin function
    takes the stream id and number of streams ("i" and "n" from the for all
    streams pragma) and returns a thread_iterator pointing to the beginning
    of the range for the stream.  The end function is similar except it
    returns a thread_iterator that is a one-past-the-end iterator for the
    stream's range. The thread_iterator is similar to an STL forward iterator.
    One additional operation provided by the thread_iterator is subtraction
    between two thread_iterators.  This returns the number of occupied entries
    in the table between the two iterators.  The below example is a function
    that takes a hash table and assigns contiguous ids (starting with 0) to
    each of the elements in the hash table.

      template<typename HT>
      inline
      void assign_contiguous_ids(HT& xht)
      {
        typename HT::size_type stream_id = 0;
        typename HT::size_type num_streams = 1;

        typename HT::size_type starting_value = 0;

        #pragma mta for all streams stream_id of num_streams
        {
          typename HT::thread_iterator begin_iter =
              xht.begin(stream_id, num_streams);

          typename HT::thread_iterator end_iter =
              xht.end(stream_id, num_streams);

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

    TODO: Add resize() and erase().  The resize() function will allow the user
          to select new table and region sizes.
*/
template <typename K, typename T, typename HF = default_hash_func<K>,
          typename EQF = default_eqfcn<K> >
class hachar {
private:
  class region;

public:
  class hash_node;

  typedef K key_type;
  typedef T mapped_type;
  typedef hash_node value_type;
  typedef EQF key_compare;
  typedef hash_size_type size_type;
  class thread_iterator;

  /// \brief The normal constructor.
  /// \param requested_table_size The size of the table.
  /// \param requested_region_size The size of each region.
  ///
  /// The requested_table_size can be any number, but it will be increased to
  /// a power of 2. The requested_region_size, if not given, will default to
  /// the requested_table_size.  The table is set up along with a first
  /// "overflow" region.
  inline hachar(size_type requested_table_size,
                size_type requested_region_size = 0);

  /// \brief The copy constructor, which performs a deep copy.
  /// \param h The hachar to be copied.
  inline hachar(const hachar& h);

  /// The destructor.
  inline ~hachar();

  /// \brief The assignment operator, which performs a deep copy.
  /// \param rhs The hachar to be copied.
  inline hachar& operator=(const hachar& rhs);

  /// \brief Returns the number of elements in the hash table.
  /// \return Number of elements in the hash table.
  ///
  /// This is a linear cost operation.
  inline size_type size() const;

  /// \brief Returns the beginning iterator for this stream.
  /// \return The beginning iterator for this stream.
  inline thread_iterator begin(size_type stream_id, size_type num_streams);

  /// Returns the ending iterator for this stream.
  /// \return The ending iterator for this stream.
  inline thread_iterator end(size_type stream_id, size_type num_streams);

  size_type get_table_size() const { return table_size; }

  /// \brief Returns true if there are no elements in the hash table.
  /// \return True if the table is empty.  Otherwise, false.
  ///
  /// This is a linear cost operation.
  size_type empty() const { return size() == 0; }

  /// \brief Returns a reference to the value associated with key.  If the key
  ///        isn't already in the table, the key is inserted, and the value is
  ///        initialized with its default constructor.  A reference to the
  ///        value of the newly created entry is returned.
  /// \param key The search key.
  /// \return Reference to the value associated with the key.
  ///
  /// This function is slightly less scalable and slightly slower than the
  /// other element access functions.  Occasionally, using it can cause other
  /// functions to have degraded performance and scaling.
  inline mapped_type& operator[](const key_type& key);

  /// \brief Tests if key exists in the table.
  /// \param key The search key.
  /// \return True if the key exists in the table.  Otherwise, false.
  inline bool member(const key_type& key) const;

  /// \brief Returns the value associated with key.
  /// \param key The search key.
  /// \param value Return value that holds the value associated with key.
  /// \return True if the key exists in the table.  Otherwise, false.
  inline bool lookup(const key_type& key, mapped_type& value) const;

  /// \brief Inserts an element into the hash table.
  /// \param key The insertion key.
  /// \param value The insertion value.
  /// \return True if the key was inserted.  Otherwise, false.
  /// \return Pair of value_type* and bool.  The first is a pointer to the
  ///         inserted key/value pair.  If the key already existed, it points
  ///         to the existing key/value pair.  The second is true if the key
  ///         was inserted and false if it already existed.
  inline pair<value_type*, bool>
  insert(const key_type& key, const mapped_type& value);

  /// \brief Inserts an element into the hash table and applies the visitor if
  ///        the key was inserted.
  /// \param key The insertion key.
  /// \param value The insertion value.
  /// \param visitor The visitor function object that is applied to the value
  ///                associated with key.
  /// \return Pair of value_type* and bool.  The first is a pointer to the
  ///         inserted key/value pair.  If the key already existed, it points
  ///         to the existing key/value pair.  The second is true if the key
  ///         was inserted and false if it already existed.
  template <typename Vis>
  inline pair<value_type*, bool>
  insert(const key_type& key, const mapped_type& value, Vis visitor);

  /// \brief Updates the value associated with key.
  /// \param key The update key.
  /// \param value The new value to associate with key.
  /// \return True if the key existed in the table.  Otherwise, false.
  inline bool update(const key_type& key, const mapped_type& value);

  /// \brief Updates the value associated with key by applying the visitor.
  /// \tparam Vis Type of the visitor function object.
  /// \param key The update key.
  /// \param value Parameter passed to the visitor function object.
  /// \param visitor The visitor function object that is applied to the value
  ///                associated with key.
  /// \return True if the key existed in the table.  Otherwise, false.
  template <typename Vis> inline
  bool update(const key_type& key, const mapped_type& value, Vis visitor);

  /// \brief Updates the value associated with key.  If key doesn't exist in
  ///        the table, insert it with a value of value.
  /// \param key The update key.
  /// \param value The new value to associate with key.
  /// \return True if the key was inserted into the table.  Otherwise, false.
  inline bool update_insert(const key_type& key, const mapped_type& value);

  /// \brief Updates the value associated with key by applying the visitor.
  ///        If key doesn't exist in the table, insert it with a value of
  ///        value.
  /// \tparam Vis Type of the visitor function object.
  /// \param key The update key.
  /// \param value Parameter passed to the visitor function object or the
  ///              value associated with key if key doesn't already exist in
  ///              the hash table.
  /// \param visitor The visitor function object that is applied to the value
  ///                associated with key.
  /// \return True if the key was inserted into the table.  Otherwise, false.
  template <typename Vis> inline
  bool update_insert(const key_type& key, const mapped_type& value,
                     Vis visitor);

  /// \brief Clears all entries from the hash table.
  ///
  /// This doesn't reduce the amount of memory allocted to the table.
  inline void clear();

  /// \brief Apply the method encapsulated by "visitor" to the key and data
  ///        of each table element, in parallel.
  /// \tparam Vis Type of the visitor function object.
  /// \param visitor A closure object containing a function to be applied to
  ///        each element of the table.
  template <typename Vis> inline void visit(Vis visitor);

  /// \brief Prints entire hash table.
  /// \param out The FILE to send the output to.
  void print(FILE* out = stdout) const;

  /// \brief Prints the number of buckets with 0, 1, etc. entries.
  /// \param out The FILE to send the output to.
  void print_stats(FILE* out = stdout) const;

  /// \brief Prints the allocated regions and the number of entries used in
  ///        each one.
  /// \param out The FILE to send the output to.
  void print_array_chunk_structure(FILE* out = stdout) const;

private:
  /// \brief Returns the index into table where key belongs.
  /// \return Index in the table where key belongs.
  ///
  /// Ands the user supplied hash function with the mask to get an index
  /// within the range 0 <= index < capacity().
  size_type hash(const key_type& key) const
  { return hash_func(key) & mask; }

  /// \brief Returns an unused node; allocates new memory, if necessary.
  /// \return A new node to put in the table.
  ///
  /// The next available node in the tail region is returned.  If the region
  /// has no more empty nodes, a new region is allocated and added to the
  /// extended table list.
  inline hash_node* allocate_link_node();

  /// \brief Performs a deep copy from rhs to the calling object.
  /// \param rhs The hachar to be copied.
  inline void deep_copy(const hachar& rhs);

  /// \brief Computes the smallest power of 2 that is greater than or equal
  ///        to the requested size.
  /// \param requested_size The requested size.
  /// \return The smallest power of 2 that is greater than or equal to the
  ///         requested size.
  size_type compute_size(size_type requested_size) const
  {
    size_type size = 4;

    if (1024 < requested_size)
    {
      size = 1024;
      if (size * size < requested_size) size *= 1024;
    }

    while (size < requested_size) size *= 2;

    return size;
  }

  size_type table_size;
  size_type region_size;
  size_type mask;
  size_type on_deck_threshold;

  hash_node* table;
  region* extended_table;
  region* tail;
  size_type* bucket_size;
  EQF compare;
  HF hash_func;
};

/// \brief This is the type of a linked list (chain) entry.
template <typename K, typename T, typename HF, typename EQF>
class hachar<K, T, HF, EQF>::hash_node {
public:
  hash_node() : first(K()), second(T()), next(NULL) {}

public:
  const K first;
  T second;

private:
  friend class hachar<K, T, HF, EQF>;
  friend class hachar<K, T, HF, EQF>::thread_iterator;

  hash_node* next;
};

/// \brief This is a forward iterator that will always be given in pairs to
///        represent a range of values in the hash table.  This iterator
///        should only be used in "for all streams i of n" sections of code.
///
/// This iterator will always point to an occupied entry in the table.  On
/// initialization, the iterator is advanced to the next occupied entry.  If
/// the end of the table is reached, the iterator will point to one past the
/// end of the table.
template <typename K, typename T, typename HF, typename EQF>
class hachar<K, T, HF, EQF>::thread_iterator {
public:
  thread_iterator() :
    table(0), bucket_size(0), data(0), index(0), table_size(0) {}

  thread_iterator(hachar& t, size_type idx) :
    table(t.table), bucket_size(t.bucket_size),
    index(idx), table_size(t.table_size)
  {
    // Advance to the first bucket that contains an item and set data to the
    // first item in that bucket.
    while (index != table_size && bucket_size[index] == 0) ++index;
    data = table + index;
  }

  thread_iterator(const thread_iterator& rhs) :
    table(rhs.table), bucket_size(rhs.bucket_size), data(rhs.data),
    index(rhs.index), table_size(rhs.table_size) {}

  thread_iterator& operator=(const thread_iterator& rhs)
  {
    table = rhs.table;
    bucket_size = rhs.bucket_size;
    data = rhs.data;
    index = rhs.index;
    table_size = rhs.table_size;
    return *this;
  }

  /// \brief Prefix increment.  Moves the iterator to the next occupied
  ///        position in the table and returns an iterator pointing to the
  ///        new position.
  thread_iterator& operator++()
  {
    do
    {
      if (data->next)
      {
        data = data->next;
      }
      else
      {
        ++index;
        data = table + index;
      }
    } while (index != table_size && bucket_size[index] == 0);

    return *this;
  }

  /// \brief Postfix increment.  Moves the iterator to the next occupied
  ///        position in the table and returns an iterator pointing to the
  ///        old position.
  thread_iterator operator++(int)
  {
    thread_iterator temp(*this);

    do
    {
      if (data->next)
      {
        data = data->next;
      }
      else
      {
        ++index;
        data = table + index;
      }
    } while (index != table_size && bucket_size[index] == 0);

    return temp;
  }

  value_type& operator*() { return *data; }
  value_type* operator->() { return data; }
  const value_type& operator*() const { return *data; }
  const value_type* operator->() const { return data; }

  bool operator==(const thread_iterator& rhs)
  { return data == rhs.data; }

  bool operator!=(const thread_iterator& rhs)
  { return data != rhs.data; }

  bool operator<(const thread_iterator& rhs)
  {
    if (index < rhs.index)
    {
      return true;
    }
    else if (index > rhs.index || data == rhs.data)
    {
      return false;
    }

    hash_node* ptr;
    for (ptr = data->next; ptr != rhs.data && ptr != NULL; ptr = ptr->next);

    return ptr != NULL;
  }

  bool operator<=(const thread_iterator& rhs)
  {
    if (index < rhs.index || data == rhs.data)
    {
      return true;
    }
    else if (index > rhs.index)
    {
      return false;
    }

    hash_node* ptr;
    for (ptr = data->next; ptr != rhs.data && ptr != NULL; ptr = ptr->next);

    return ptr != NULL;
  }

  bool operator>(const thread_iterator& rhs)
  {
    if (index > rhs.index)
    {
      return true;
    }
    else if (index < rhs.index || data == rhs.data)
    {
      return false;
    }

    hash_node* ptr;
    for (ptr = data->next; ptr != rhs.data && ptr != NULL; ptr = ptr->next);

    return ptr == NULL;
  }

  bool operator>=(const thread_iterator& rhs)
  {
    if (index > rhs.index || data == rhs.data)
    {
      return true;
    }
    else if (index < rhs.index)
    {
      return false;
    }

    hash_node* ptr;
    for (ptr = data->next; ptr != rhs.data && ptr != NULL; ptr = ptr->next);

    return ptr == NULL;
  }

  /// Returns the number of occupied table entries between the two iterators.
  size_type operator-(const thread_iterator& rhs)
  {
    if (data == rhs.data) return 0;

    // Account for rhs.
    size_type num_entries = 1;

    // Count the entries after rhs in the first bucket.
    hash_node* ptr;
    for (ptr = rhs.data->next; ptr != NULL; ptr = ptr->next) ++num_entries;

    // Count the entries in the middle buckets.
    size_type i;
    for (i = rhs.index + 1; i < index; ++i)
    {
      num_entries += bucket_size[i];
    }

    // Count the entries in the last bucket up to but not including lhs.
    for (ptr = table + i; ptr != data; ptr = ptr->next) ++num_entries;

    return num_entries;
  }

private:
  hash_node* table;
  size_type* bucket_size;
  hash_node* data;
  size_type index;
  size_type table_size;
};

/// \brief This is the type of a region in the hachar structure.
template <typename K, typename T, typename HF, typename EQF>
class hachar<K, T, HF, EQF>::region {
public:
  /// \brief The normal constructor.
  /// \param size The size of this region.
  region(size_type chunkSize) : chunk_used(0), next(NULL)
  {
    chunk = (hash_node*) malloc(chunkSize * sizeof(hash_node));
  }

  /// \brief The destructor.
  ~region() { free(chunk); }

  hash_node* chunk;
  size_type chunk_used;
  region* next;
};

template <typename K, typename T, typename HF, typename EQF>
hachar<K, T, HF, EQF>::hachar(size_type requested_table_size,
                              size_type requested_region_size)
{
  table_size = compute_size(requested_table_size);
  region_size = requested_region_size == 0 ? table_size : requested_region_size;
  mask = table_size - 1;
  on_deck_threshold = region_size / 2;

  table = (hash_node*) malloc(table_size * sizeof(hash_node));
  extended_table = new region(region_size);
  bucket_size = (size_type*) malloc(table_size * sizeof(size_type));
  tail = extended_table;

  size_type tsize = table_size;
  #pragma mta assert nodep
  for (size_type i = 0; i < tsize; ++i)
  {
    table[i].next = NULL;
    bucket_size[i] = 0;
  }
}

template <typename K, typename T, typename HF, typename EQF>
hachar<K, T, HF, EQF>::hachar(const hachar& h)
{
  deep_copy(h);
}

template <typename K, typename T, typename HF, typename EQF>
hachar<K, T, HF, EQF>::~hachar()
{
  free(bucket_size);
  free(table);

  while (extended_table != NULL)
  {
    region* ptr = extended_table->next;
    delete extended_table;
    extended_table = ptr;
  }
};

template <typename K, typename T, typename HF, typename EQF>
hachar<K, T, HF, EQF>&
hachar<K, T, HF, EQF>::operator=(const hachar& rhs)
{
  free(bucket_size);
  free(table);

  while (extended_table != NULL)
  {
    region* ptr = extended_table->next;
    delete extended_table;
    extended_table = ptr;
  }

  deep_copy(rhs);

  return *this;
}

template <typename K, typename T, typename HF, typename EQF>
typename hachar<K, T, HF, EQF>::thread_iterator
hachar<K, T, HF, EQF>::begin(size_type stream_id, size_type num_streams)
{
  size_type index = begin_block_range(table_size, stream_id, num_streams);
  return thread_iterator(*this, index);
}

template <typename K, typename T, typename HF, typename EQF>
typename hachar<K, T, HF, EQF>::thread_iterator
hachar<K, T, HF, EQF>::end(size_type stream_id, size_type num_streams)
{
  size_type index = end_block_range(table_size, stream_id, num_streams);
  return thread_iterator(*this, index);
}

template <typename K, typename T, typename HF, typename EQF>
typename hachar<K, T, HF, EQF>::size_type
hachar<K, T, HF, EQF>::size() const
{
  size_type result = 0;
  size_type tsize = table_size;

  #pragma mta assert nodep
  for (size_type i = 0; i < tsize; ++i) result += bucket_size[i];

  return result;
}

template <typename K, typename T, typename HF, typename EQF>
T&
hachar<K, T, HF, EQF>::operator[](const K& key)
{
  size_type index = hash(key);

  size_type& counter = bucket_size[index];
  hash_node* ptr = table + index;

  // Insert in the main table if that entry is empty.
  if (counter == 0)
  {
    long probed = mt_readfe(counter);

    if (probed == 0)
    {
      const_cast<K&>(ptr->first) = key;
      ptr->second = T();
      mt_write(counter, 1);
      return ptr->second;
    }

    mt_write(counter, probed);
  }

  // Search the bucket's linked list.  If the key is found, return its value.
  // If the end of the chain is reached, the key doesn't exist.  Add a new
  // entry to the end of the linked list, insert the key, and return the
  // default constructor initialized value.
  while (compare(ptr->first, key) == false)
  {
    #pragma mta expect true
    if (ptr->next == NULL)
    {
      // The key isn't in the bucket.  Add a new entry to the bucket.
      hash_node* probed = mt_readfe(ptr->next);

      if (probed == NULL)
      {
        hash_node* q = allocate_link_node();
        q->next = NULL;
        const_cast<K&>(q->first) = key;
        q->second = T();
        mt_write(ptr->next, q);
        mt_incr(counter, 1);
        return q->second;
      }

      mt_write(ptr->next, probed);
    }
    else
    {
      // Search the next entry in the bucket.
      ptr = ptr->next;
    }
  }

  // The key was found.  Return its value.
  return ptr->second;
}

template <typename K, typename T, typename HF, typename EQF>
bool
hachar<K, T, HF, EQF>::member(const K& key) const
{
  size_type index = hash(key);

  if (bucket_size[index] == 0) return false;

  hash_node* ptr = table + index;
  do
  {
    if (compare(ptr->first, key)) return true;

    ptr = ptr->next;
  } while (ptr != NULL);

  return false;
}

template <typename K, typename T, typename HF, typename EQF>
bool
hachar<K, T, HF, EQF>::lookup(const K& key, T& value) const
{
  size_type index = hash(key);

  if (bucket_size[index] == 0) return false;

  hash_node* ptr = table + index;
  do
  {
    if (compare(ptr->first, key))
    {
      value = ptr->second;
      return true;
    }

    ptr = ptr->next;
  } while (ptr != NULL);

  return false;
}

template <typename K, typename T, typename HF, typename EQF>
pair<typename hachar<K, T, HF, EQF>::value_type*, bool>
hachar<K, T, HF, EQF>::insert(const K& key, const T& value)
{
  size_type index = hash(key);

  size_type& counter = bucket_size[index];
  hash_node* ptr = table + index;

  // Insert in the main table if that entry is empty.
  if (counter == 0)
  {
    long probed = mt_readfe(counter);

    if (probed == 0)
    {
      const_cast<K&>(ptr->first) = key;
      ptr->second = value;
      mt_write(counter, 1);
      return pair<value_type*, bool>(ptr, true);
    }

    mt_write(counter, probed);
  }

  // Search the bucket's linked list.  If the key is found, just return false.
  // If the end of the chain is reached, the key doesn't exist.  Add a new
  // entry to the end of the linked list and insert the key.
  while (compare(ptr->first, key) == false)
  {
    #pragma mta expect true
    if (ptr->next == NULL)
    {
      // The key isn't in the bucket.  Add a new entry to the bucket.
      hash_node* probed = mt_readfe(ptr->next);

      if (probed == NULL)
      {
        hash_node* q = allocate_link_node();
        q->next = NULL;
        const_cast<K&>(q->first) = key;
        q->second = value;
        mt_write(ptr->next, q);
        mt_incr(counter, 1);
        return pair<value_type*, bool>(q, true);
      }

      mt_write(ptr->next, probed);
    }
    else
    {
      // Search the next entry in the bucket.
      ptr = ptr->next;
    }
  }

  // The key was found.  Return that the key wasn't inserted.
  return pair<value_type*, bool>(ptr, false);
}

template <typename K, typename T, typename HF, typename EQF>
template <typename Vis>
pair<typename hachar<K, T, HF, EQF>::value_type*, bool>
hachar<K, T, HF, EQF>::insert(const K& key, const T& value, Vis visitor)
{
  size_type index = hash(key);

  size_type& counter = bucket_size[index];
  hash_node* ptr = table + index;

  // Insert in the main table if that entry is empty, and apply the visitor.
  if (counter == 0)
  {
    long probed = mt_readfe(counter);

    if (probed == 0)
    {
      const_cast<K&>(ptr->first) = key;
      ptr->second = value;
      visitor(ptr->second);
      mt_write(counter, 1);
      return pair<value_type*, bool>(ptr, true);
    }

    mt_write(counter, probed);
  }

  // Search the bucket's linked list.  If the key is found, just return false.
  // If the end of the chain is reached, the key doesn't exist.  Add a new
  // entry to the end of the linked list, insert the key, and apply the
  // visitor.
  while (compare(ptr->first, key) == false)
  {
    #pragma mta expect true
    if (ptr->next == NULL)
    {
      // The key isn't in the bucket.  Add a new entry to the bucket.
      hash_node* probed = mt_readfe(ptr->next);

      if (probed == NULL)
      {
        hash_node* q = allocate_link_node();
        q->next = NULL;
        const_cast<K&>(q->first) = key;
        q->second = value;
        visitor(q->second);
        mt_write(ptr->next, q);
        mt_incr(counter, 1);
        return pair<value_type*, bool>(q, true);
      }

      mt_write(ptr->next, probed);
    }
    else
    {
      // Search the next entry in the bucket.
      ptr = ptr->next;
    }
  }

  // The key was found.  Return that the key wasn't inserted.
  return pair<value_type*, bool>(ptr, false);
}

template <typename K, typename T, typename HF, typename EQF>
bool
hachar<K, T, HF, EQF>::update(const K& key, const T& value)
{
  size_type index = hash(key);

  if (bucket_size[index] == 0) return false;

  hash_node* ptr = table + index;
  do
  {
    if (compare(ptr->first, key))
    {
      ptr->second = value;
      return true;
    }

    ptr = ptr->next;
  } while (ptr != NULL);

  return false;
}

template <typename K, typename T, typename HF, typename EQF>
template <typename Vis>
bool
hachar<K, T, HF, EQF>::update(const K& key, const T& value, Vis visitor)
{
  size_type index = hash(key);

  if (bucket_size[index] == 0) return false;

  hash_node* ptr = table + index;
  do
  {
    if (compare(ptr->first, key))
    {
      visitor(ptr->second, value);
      return true;
    }

    ptr = ptr->next;
  } while (ptr != NULL);

  return false;
}

template <typename K, typename T, typename HF, typename EQF>
bool
hachar<K, T, HF, EQF>::update_insert(const K& key, const T& value)
{
  size_type index = hash(key);

  size_type& counter = bucket_size[index];
  hash_node* ptr = table + index;

  // Insert in the main table if that entry is empty.
  if (counter == 0)
  {
    long probed = mt_readfe(counter);

    if (probed == 0)
    {
      const_cast<K&>(ptr->first) = key;
      ptr->second = value;
      mt_write(counter, 1);
      return true;
    }

    mt_write(counter, probed);
  }

  // Search the bucket's linked list.  If the key is found, update the key.
  // If the end of the chain is reached, the key doesn't exist.  Add a new
  // entry to the end of the linked list and insert the key.
  while (compare(ptr->first, key) == false)
  {
    #pragma mta expect false
    if (ptr->next == NULL)
    {
      // The key isn't in the bucket.  Add a new entry to the bucket.
      hash_node* probed = mt_readfe(ptr->next);

      if (probed == NULL)
      {
        hash_node* q = allocate_link_node();
        q->next = NULL;
        const_cast<K&>(q->first) = key;
        q->second = value;
        mt_write(ptr->next, q);
        mt_incr(counter, 1);
        return true;
      }

      mt_write(ptr->next, probed);
    }
    else
    {
      // Search the next entry in the bucket.
      ptr = ptr->next;
    }
  }

  // The key was found.  Update its value.
  ptr->second = value;
  return false;
}

template <typename K, typename T, typename HF, typename EQF>
template <typename Vis>
bool
hachar<K, T, HF, EQF>::update_insert(const K& key, const T& value, Vis visitor)
{
  size_type index = hash(key);

  size_type& counter = bucket_size[index];
  hash_node* ptr = table + index;

  // Insert in the main table if that entry is empty.
  if (counter == 0)
  {
    long probed = mt_readfe(counter);

    if (probed == 0)
    {
      const_cast<K&>(ptr->first) = key;
      ptr->second = value;
      mt_write(counter, 1);
      return true;
    }

    mt_write(counter, probed);
  }

  // Search the bucket's linked list.  If the key is found, apply the visitor
  // to its value.  If the end of the chain is reached, the key doesn't exist.
  // Add a new entry to the end of the linked list and insert the key.
  while (compare(ptr->first, key) == false)
  {
    #pragma mta expect false
    if (ptr->next == NULL)
    {
      // The key isn't in the bucket.  Add a new entry to the bucket.
      hash_node* probed = mt_readfe(ptr->next);

      if (probed == NULL)
      {
        hash_node* q = allocate_link_node();
        q->next = NULL;
        const_cast<K&>(q->first) = key;
        q->second = value;
        mt_write(ptr->next, q);
        mt_incr(counter, 1);
        return true;
      }

      mt_write(ptr->next, probed);
    }
    else
    {
      // Search the next entry in the bucket.
      ptr = ptr->next;
    }
  }

  // The key was found.  Apply the visitor to its value.
  visitor(ptr->second, value);
  return false;
}

template <typename K, typename T, typename HF, typename EQF>
void
hachar<K, T, HF, EQF>::clear()
{
  // Set the chunk_used of all regions except the main table to 0.
  region* region_ptr = extended_table;
  while (region_ptr != NULL)
  {
    region_ptr->chunk_used = 0;
    region_ptr = region_ptr->next;
  }

  // Set the values in bucket_size to 0.
  memset(bucket_size, 0, table_size * sizeof(size_type));

  size_type tsize = table_size;

  // Follow the linked lists in the buckets setting all pointers to NULL.
  #pragma mta assert nodep
  for (size_type i = 0; i < tsize; ++i)
  {
    hash_node* ptr = table + i;

    while (ptr->next != NULL)
    {
      hash_node* tmp_ptr = ptr->next;
      ptr->next = NULL;
      ptr = tmp_ptr;
    }
  }
}

template <typename K, typename T, typename HF, typename EQF>
template <typename Vis>
void
hachar<K, T, HF, EQF>::visit(Vis visitor)
{
  size_type tsize = table_size;

  // Follow the linked lists in the buckets applying the visitor to each
  // entry.
  #pragma mta assert parallel
  for (size_type i = 0; i < tsize; ++i)
  {
    hash_node* ptr = table + i;

    while (ptr != NULL)
    {
      visitor(ptr);
      ptr = ptr->next;
    }
  }
}

template <typename K, typename T, typename HF, typename EQF>
typename hachar<K, T, HF, EQF>::hash_node*
hachar<K, T, HF, EQF>::allocate_link_node()
{
  region* ptr = mt_readxx(tail);
  size_type index = mt_incr(ptr->chunk_used, 1);

  #pragma mta expect false
  if (index == on_deck_threshold && mt_readxx(ptr->next) == NULL)
  {
    mt_writexf(ptr->next, new region(region_size));
  }

  while (index >= region_size)
  {
    #pragma mta expect false
    if (index == region_size)
    {
      if (mt_readxx(ptr->next) == NULL) continue;

      mt_writexf(tail, mt_readxx(ptr->next));
    }

    ptr = mt_readxx(tail);
    index = mt_incr(ptr->chunk_used, 1);
  }

  return &(ptr->chunk[index]);
}

template <typename K, typename T, typename HF, typename EQF>
void
hachar<K, T, HF, EQF>::deep_copy(const hachar& rhs)
{
  table_size = rhs.table_size;
  region_size = rhs.region_size;
  mask = rhs.mask;
  on_deck_threshold = rhs.on_deck_threshold;

  table = (hash_node*) malloc(table_size * sizeof(hash_node));
  extended_table = new region(region_size);
  bucket_size = (size_type*) malloc(table_size * sizeof(size_type));
  tail = extended_table;

  size_type tsize = table_size;

  // Copy the bucket size and intialize the next pointer of the table nodes.
  #pragma mta assert nodep
  for (size_type i = 0; i < tsize; ++i) bucket_size[i] = rhs.bucket_size[i];

  // Copy the entries in the table.
  #pragma mta assert nodep
  #pragma mta block dynamic schedule
  for (size_type i = 0; i < tsize; ++i)
  {
    if (bucket_size[i] > 0)
    {
      const_cast<K&>(table[i].first) = rhs.table[i].first;
      table[i].second = rhs.table[i].second;
    }
  }

  // Copy the entries in the extended table by following the linked lists in
  // the buckets.
  #pragma mta assert parallel
  #pragma mta block dynamic schedule
  for (size_type i = 0; i < tsize; ++i)
  {
    if (bucket_size[i] > 1)
    {
      hash_node* ptr = table + i;
      hash_node* rhs_ptr = rhs.table[i].next;

      while (rhs_ptr != NULL)
      {
        ptr->next = allocate_link_node();
        ptr = ptr->next;
        const_cast<K&>(ptr->first) = rhs_ptr->first;
        ptr->second = rhs_ptr->second;
        rhs_ptr = rhs_ptr->next;
      }

      ptr->next = NULL;
    }
    else
    {
      table[i].next = NULL;
    }
  }
}

template <typename K, typename T, typename HF, typename EQF>
void hachar<K, T, HF, EQF>::print(FILE* out) const
{
  fprintf(out, "HashTable Structure:\n");
  region* region_ptr = extended_table;

  while (region_ptr != NULL)
  {
    fprintf(out, " ---> buffer(%ps) contains %lu items\n", region_ptr,
            region_ptr->chunk_used);
    region_ptr = region_ptr->next;
  }

  fprintf(out, " -> Bucket Sizes:\n");

  for (size_type i = 0; i < table_size; ++i)
  {
    fprintf(out, "%5lu: %2lu", i, bucket_size[i]);

    if (bucket_size[i] > 0)
    {
      hash_node* ptr = table + i;

      while (ptr != NULL)
      {
        std::ostringstream ostr1;
        std::ostringstream ostr2;
        ostr1 << ptr->first;
        ostr2 << ptr->second;
        fprintf(out, "  (%s => %s)", ostr1.str().c_str(), ostr2.str().c_str());
        ptr = ptr->next;
      }
    }

    fprintf(out, "\n");
  }
}

template <typename K, typename T, typename HF, typename EQF>
void hachar<K, T, HF, EQF>::print_stats(FILE* out) const
{
  size_type tsize = table_size;
  size_type num_entries = size();

  fprintf(out, "Hash table slots: %lu, number of entries: %lu\n",
          tsize, num_entries);

  size_type maxsz = 0;
  for (size_type i = 0; i < tsize; ++i)
  {
    if (bucket_size[i] > maxsz) maxsz = bucket_size[i];
  }

  size_type histsz = maxsz + 1;
  size_type* hist = new size_type[histsz];

  for (size_type i = 0; i < histsz; ++i) hist[i] = 0;

  for (size_type i = 0; i < tsize; ++i) ++hist[bucket_size[i]];

  for (size_type i = 0; i < histsz; ++i)
  {
    fprintf(out, "Number of %lu: %lu\n", i, hist[i]);
  }

  delete [] hist;
}

template <typename K, typename T, typename HF, typename EQF>
void hachar<K, T, HF, EQF>::print_array_chunk_structure(FILE* out) const
{
  fprintf(out, "ArrayChunk structure:\n");
  if (tail->next != NULL) fprintf(out, "  -> on deck: %p\n", tail->next);

  region* region_ptr = extended_table;
  while (region_ptr != NULL)
  {
    fprintf(out, "  -> buffer: %p, chunk_used=%lu\n", region_ptr,
            region_ptr->chunk_used);
    region_ptr = region_ptr->next;
  }
}

}

#endif
