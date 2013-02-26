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
/*! \file xmt_hash_table.hpp

    \brief A thread-safe key-based dictionary abstraction implemented using
           open address hashing.

    \author Greg Mackey (gemacke@sandia.gov)
    \author Eric Goodman (elgoodm@sandia.gov)
    \author Jon Berry (jberry@sandia.gov)

    \date 4/5/2008
*/
/****************************************************************************/

#ifndef MTGL_XMT_HASH_TABLE_HPP
#define MTGL_XMT_HASH_TABLE_HPP

#include <limits>
#include <iostream>

#include <mtgl/util.hpp>
#include <mtgl/partitioning.hpp>
#include <mtgl/hash_defs.hpp>
#include <mtgl/mmap_traits.hpp>

#define XMT_HT_EMPTY   0
#define XMT_HT_DELETED 1
#define XMT_HT_OCCUPIED 2

namespace mtgl {

/*! \brief A map data structure implemented using hashing with open
           addressing.

    \author Greg Mackey (gemacke@sandia.gov)
    \author Eric Goodman (elgoodm@sandia.gov)
    \author Jon Berry (jberry@sandia.gov)

    \tparam K Key type.
    \tparam T Value type.
    \tparam HF Function object type defining the hash function.
    \tparam EQF Function object type defining the equality function for keys.

    The hash table uses fine-grained locking at the entry level, and it uses
    the double locking method described in the paper "Hashing Strategies for
    the Cray XMT" by Goodman, et. al.  It also avoids locking whenever
    possible.  Locking is only used when modifying a key.  All of these
    methods increase the parallelism of the hash table.

    The hash table enforces unique keys, and it allows deletes.  Note that
    keys cannot be modified.  Once they are inserted into the hash, they can
    only be deleted.  The values, however, are modifiable via the update or
    the visit functions.  Because performing an int_fetch_add during insert on
    an element counter causes a hotspot, the table doesn't keep track of
    the number of elements currently in the table.  Instead, a linear cost
    size() function is provided.

    The [] operator is slightly less scalable and slightly slower than the
    other element access functions (insert, update, etc.).  Those functions
    should be preferred over the [] operator when possible.

    The hash table performs deletions by marking the entries as deleted.  The
    normal insert functions ignore these deleted entries when searching for an
    open spot in the table.  This method of deletion is known as tombstoning.
    Once a position in the table is deleted, the position is unavailable for
    use by later insertions.  Thus, erasing elements shrinks the size of the
    table.  This is done because reusing the elements makes the inserts
    more expensive in instruction count, loads, and stores, although the time
    performance was the same in my tests for inserts into an empty table.  If
    space is more of a concern than performance, insert functions that reuse
    the deleted positions are provided.  They are the insert functions whose
    names end in "_reuse".

    There are two main classes of functions: those designed to be called in a
    parallel context and those designed to be called in a serial context.  The
    parallel context functions can be further divided into functions that
    modify the keys and functions that don't.  Here is the division of
    functions:

      - Parallel Context Functions
         - Modify the Key
            - operator[]()
            - insert()
            - insert_reuse()
            - erase()
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
         - resize()
         - num_deleted_entries()
         - clear_deleted_entries()
         - capacity()
         - empty()
         - swap()
         - clear()
         - visit()
         - print()

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
    and not around modifying the value.  As an example, consider the case
    where one thread is updating an entry and another thread is deleting the
    same entry.  The expected outcome is that the entry is deleted.  Let the
    update thread find the entry.  The delete thread then completes its
    delete.  The update thread then updates the value associated with the
    deleted entry.  This is okay because the entry is still deleted.

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

      template <typename HT>
      class table_visitor {
      public:
        void operator()(typename HT::value_type* i) { i->second += i->first; }
      };

    Now an example of using the visit function.

      xmt_hash_table<int, int> ht;
      <Items inserted into the hash here.>
      table_visitor<xmt_hash_table<int, int> > tv;
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

    The hash table offers static resizing via the resize function.  Note that
    this function is the only way to resize the table.  Also, it is not
    thread-safe in the sense that other operations on the hash table cannot be
    called at the same time as the resize function.  The user must ensure that
    the resize function is not interleaved with other hash table functions.
    (i.e. The user shouldn't put the resize function inside a loop.)  Making
    the resize function completely thread-safe would require implementing a
    reader-writer locking scheme which could significantly slow down the hash
    table performance.
*/
template <typename K, typename T, typename HF = default_hash_func<K>,
          typename EQF = default_eqfcn<K> >
class xmt_hash_table {
public:
  typedef K key_type;
  typedef T mapped_type;
  typedef pair<const K, T> value_type;
  typedef EQF key_compare;
  typedef hash_size_type size_type;
  class thread_iterator;

  /// \brief The normal constructor.
  /// \param size The suggested initial allocation size.
  ///
  /// The allocation size must always be a power of two, so the allocation
  /// size is set to the smallest power of two that is greater than or equal
  /// to size.
  inline xmt_hash_table(size_type size = 1024);

  /// \brief The copy constructor, which performs a deep copy.
  /// \param h The xmt_hash_table to be copied.
  inline xmt_hash_table(const xmt_hash_table& h);

  /// The destructor.
  inline ~xmt_hash_table();

  /// \brief The assignment operator, which performs a deep copy.
  /// \param rhs The xmt_hash_table to be copied.
  inline xmt_hash_table& operator=(const xmt_hash_table& rhs);

  /// Returns the beginning iterator for this stream.
  inline thread_iterator begin(size_type stream_id, size_type num_streams);

  /// Returns the ending iterator for this stream.
  inline thread_iterator end(size_type stream_id, size_type num_streams);

  /// \brief Returns the number of elements in the hash table.
  /// \return Number of elements in the hash table.
  ///
  /// This is a linear cost operation.
  inline size_type size() const;

  /// \brief Resizes the hash table.  Any existing entries are rehased to the
  ///        new table.
  /// \param new_size The suggested new size for the hash table.
  ///
  /// The allocation size must always be a power of two, and the new hash
  /// table size must be at least the number of entries currently in the
  /// table.  Thus, the allocation size is set to the smallest power of two
  /// that is greater than or equal to the larger of the current number of
  /// entries and new_size.
  ///
  /// Note that calling resize on an mmapped xmt_hash_table does not change
  /// the data located in the mmapped region of memory.  It creates new memory
  /// and rehashes the entries in the mmapped region to that new memory.
  inline void resize(size_type new_size);

  /// \brief Returns the number of deleted entries in the hash table.
  /// \return Number of deleted entries in the hash table.
  ///
  /// This is a linear cost operation.  When the number of the deleted
  /// in the table is high, the cost of inserting into the table becomes
  /// expensive.  This lets you test if you have too many deleted entries
  /// and need to clear the deleted entries to improve the performance of
  /// the table.
  inline size_type num_deleted_entries() const;

  /// \brief Clears the deleted entries from the hash table to improve
  ///        performance.
  ///
  /// This function is relatively expensive because cleaning out the deleted
  /// entries requires rehashing all the exisitng elements into new memory.
  /// After calling this function, all the entries in the new table are either
  /// occupied or empty.
  inline void clear_deleted_entries() { resize(table_size); }

  /// \brief Returns the allocation size of the table.
  /// \return Allocation size of table.
  size_type capacity() const { return table_size; }

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
  /// \return Pair of value_type* and bool.  The first is a pointer to the
  ///         inserted key/value pair.  If the key already existed, it points
  ///         to the existing key/value pair.  If the table was full, it is
  ///         set to NULL.  The second is true if the key was inserted and
  ///         false if it already existed.
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
  ///         to the existing key/value pair.  If the table was full, it is
  ///         set to NULL.  The second is true if the key was inserted and
  ///         false if it already existed.
  template <typename Vis>
  inline pair<value_type*, bool>
  insert(const key_type& key, const mapped_type& value, Vis visitor);

  /// \brief Inserts an element into the hash table reusing deleted entries.
  /// \param key The insertion key.
  /// \param value The insertion value.
  /// \return Pair of value_type* and bool.  The first is a pointer to the
  ///         inserted key/value pair.  If the key already existed, it points
  ///         to the existing key/value pair.  If the table was full, it is
  ///         set to NULL.  The second is true if the key was inserted and
  ///         false if it already existed.
  ///
  /// This version of insert is slightly more expensive than insert().
  inline pair<value_type*, bool>
  insert_reuse(const key_type& key, const mapped_type& value);

  /// \brief Inserts an element into the hash table reusing deleted entries.
  ///        Applies the visitor if the key was inserted.
  /// \param key The insertion key.
  /// \param value The insertion value.
  /// \param visitor The visitor function object that is applied to the value
  ///                associated with key.
  /// \return Pair of value_type* and bool.  The first is a pointer to the
  ///         inserted key/value pair.  If the key already existed, it points
  ///         to the existing key/value pair.  If the table was full, it is
  ///         set to NULL.  The second is true if the key was inserted and
  ///         false if it already existed.
  ///
  /// This version of insert is slightly more expensive than insert().
  template <typename Vis>
  inline pair<value_type*, bool>
  insert_reuse(const key_type& key, const mapped_type& value, Vis visitor);

  /// \brief Erases key from the table.
  /// \param key The delete key.
  /// \return True if the key existed in the table and was erased.
  ///         Otherwise, false.
  inline bool erase(const key_type& key);

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

  /// Swaps the contents of the two hash tables.
  inline void swap(xmt_hash_table& rhs);

  /// Clears all entries from the hash table.
  inline void clear();

  /// \brief Apply the method encapsulated by "visitor" to the key and data
  ///        of each table element, in parallel.
  /// \tparam Vis Type of the visitor function object.
  /// \param visitor A closure object containing a function to be applied to
  ///        each element of the table.
  template <typename Vis> inline void visit(Vis visitor);

  /// Prints the key-value pairs in the hash table.
  void print() const;

  /// \brief Returns the amount of contiguous memory needed to mmap this
  ///        hash table.
  unsigned long get_mmap_size();

  /// \brief Writes the data for the hash table to a previously mmapped
  ///        region of memory.
  /// \param mapped_mem A pointer to the previously mmapped memory.
  void write_mmap(void* mapped_mem);

  /// \brief Clears the hash table and initializes it by pointing it to data
  ///        from another hash table that has beeen previously mmapped.
  /// \param mapped_mem A pointer to the previously mmapped memory.
  void read_mmap(void* mapped_mem);

private:
  /// \brief Performs a deep copy from rhs to the calling object.
  /// \param rhs The xmt_hash_table to be copied.
  inline void deep_copy(const xmt_hash_table& rhs);

  /// \brief Deletes all the memory allocated by the hash table.
  inline void clear_mem();

  /// \brief Returns the index into table where key belongs.
  ///
  /// Ands the user supplied hash function with the mask to get an index
  /// within the range 0 <= index < capacity().
  size_type hash(const key_type& key) const { return hash_func(key) & mask; }

private:
  size_type mask;
  size_type table_size;
  value_type* table;
  long* occupied;
  bool initialized_by_mmap;
  EQF compare;
  HF hash_func;
};

/***/

/// \brief This is a forward iterator that will always be given in pairs to
///        represent a range of values in the hash table.  This iterator
///        should only be used in "for all streams i of n" sections of code.
///
/// This iterator will always point to an occupied entry in the table.  On
/// initialization, the iterator is advanced to the next occupied entry.  If
/// the end of the table is reached, the iterator will point to one past the
/// end of the table.
template <typename K, typename T, typename HF, typename EQF>
class xmt_hash_table<K, T, HF, EQF>::thread_iterator {
public:
  thread_iterator() : occupied(0), table(0), table_size(0), index(0) {}

  thread_iterator(xmt_hash_table& t, size_type idx) :
    occupied(t.occupied), table(t.table),
    table_size(t.table_size), index(idx)
  {
    while (index != table_size && occupied[index] != XMT_HT_OCCUPIED) ++index;
  }

  thread_iterator(const thread_iterator& rhs) :
    occupied(rhs.occupied), table(rhs.table),
    table_size(rhs.table_size), index(rhs.index) {}

  thread_iterator& operator=(const thread_iterator& rhs)
  {
    occupied = rhs.occupied;
    table = rhs.table;
    table_size = rhs.table_size;
    index = rhs.index;
    return *this;
  }

  /// \brief Prefix increment.  Moves the iterator to the next occupied
  ///        position in the table and returns an iterator pointing to the
  ///        new position.
  thread_iterator& operator++()
  {
    do
    {
      ++index;
    } while (index != table_size && occupied[index] != XMT_HT_OCCUPIED);

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
      ++index;
    } while (index != table_size && occupied[index] != XMT_HT_OCCUPIED);

    return temp;
  }

  value_type& operator*() { return table[index]; }
  value_type* operator->() { return &(table[index]); }
  const value_type& operator*() const { return table[index]; }
  const value_type* operator->() const { return &(table[index]); }

  bool operator==(const thread_iterator& rhs) { return index == rhs.index; }
  bool operator!=(const thread_iterator& rhs) { return index != rhs.index; }
  bool operator<(const thread_iterator& rhs) { return index < rhs.index; }
  bool operator<=(const thread_iterator& rhs) { return index <= rhs.index; }
  bool operator>(const thread_iterator& rhs) { return index > rhs.index; }
  bool operator>=(const thread_iterator& rhs) { return index >= rhs.index; }

  /// Returns the number of occupied table entries between the two iterators.
  size_type operator-(const thread_iterator& rhs)
  {
    size_type num_entries = 0;

    for (size_type i = rhs.index; i < index; ++i)
    {
      if (occupied[i] == XMT_HT_OCCUPIED) ++num_entries;
    }

    return num_entries;
  }

private:
  long* occupied;
  value_type* table;
  size_type table_size;
  size_type index;
};

/***/

template <typename K, typename T, typename HF, typename EQF>
xmt_hash_table<K, T, HF, EQF>::xmt_hash_table(size_type size)
{
  initialized_by_mmap = false;

  size_type value = 16;
  if (size < 16)  size = 16;
  while (value < size)  { value *= 2; }
  table_size = value;
  mask = value - 1;

  // Allocate memory for table and occupied vector, and initialize occupied
  // vector to indicate all entries are empty.
  #pragma mta assert par_newdelete
  table = new value_type[table_size];
  occupied = (long*) malloc(table_size * sizeof(long));

  #pragma mta assert parallel
  #pragma mta block schedule
  for (size_type i = 0; i < table_size; ++i) occupied[i] = XMT_HT_EMPTY;
}

/***/

template <typename K, typename T, typename HF, typename EQF>
xmt_hash_table<K, T, HF, EQF>::xmt_hash_table(const xmt_hash_table& h) :
    mask(h.mask), table_size(h.table_size), initialized_by_mmap(false)
{
  deep_copy(h);
}

/***/

template <typename K, typename T, typename HF, typename EQF>
xmt_hash_table<K, T, HF, EQF>::~xmt_hash_table()
{
  clear_mem();
}

/***/

template <typename K, typename T, typename HF, typename EQF>
xmt_hash_table<K, T, HF, EQF>&
xmt_hash_table<K, T, HF, EQF>::operator=(const xmt_hash_table& rhs)
{
  clear_mem();

  mask = rhs.mask;
  table_size = rhs.table_size;

  deep_copy(rhs);

  return *this;
}

/***/

template <typename K, typename T, typename HF, typename EQF>
typename xmt_hash_table<K, T, HF, EQF>::thread_iterator
xmt_hash_table<K, T, HF, EQF>::begin(size_type stream_id, size_type num_streams)
{
  size_type index = begin_block_range(table_size, stream_id, num_streams);
  return thread_iterator(*this, index);
}

/***/

template <typename K, typename T, typename HF, typename EQF>
typename xmt_hash_table<K, T, HF, EQF>::thread_iterator
xmt_hash_table<K, T, HF, EQF>::end(size_type stream_id, size_type num_streams)
{
  size_type index = end_block_range(table_size, stream_id, num_streams);
  return thread_iterator(*this, index);
}

/***/

template <typename K, typename T, typename HF, typename EQF>
typename xmt_hash_table<K, T, HF, EQF>::size_type
xmt_hash_table<K, T, HF, EQF>::size() const
{
  size_type num_elements = 0;
  for (size_type i = 0; i < table_size; ++i)
  {
    num_elements += (occupied[i] == XMT_HT_OCCUPIED);
  }
  return num_elements;
}

/***/

template <typename K, typename T, typename HF, typename EQF>
void
xmt_hash_table<K, T, HF, EQF>::resize(size_type new_size)
{
  // Get copy of old values.
  value_type* old_table = table;
  long* old_occupied = occupied;
  size_type old_table_size = table_size;

  size_type num_elements = size();

  // Find the new size of the table.
  size_type value = 16;
  if (new_size < 16)  new_size = 16;
  while (value < new_size || value < num_elements)  { value *= 2; }
  table_size = value;
  mask = value - 1;

  // Allocate memory for table and occupied vector, and initialize occupied
  // vector to indicate all entries are empty.
  #pragma mta assert par_newdelete
  table = new value_type[table_size];
  occupied = (long*) malloc(table_size * sizeof(long));

  #pragma mta assert parallel
  #pragma mta block schedule
  for (size_type i = 0; i < table_size; ++i) occupied[i] = XMT_HT_EMPTY;

  // If the table isn't empty, copy the existing entries.
  if (num_elements > 0)
  {
    // Insert all the existing entries into the new table.
    #pragma mta assert parallel
    #pragma mta block schedule
    for (size_type i = 0; i < old_table_size; ++i)
    {
      if (old_occupied[i] == XMT_HT_OCCUPIED)
      {
        insert(old_table[i].first, old_table[i].second);
      }
    }
  }

  // Delete memory for the old table and occupied vector.
  #pragma mta assert par_newdelete
  delete [] old_table;
  free(old_occupied);
}

/***/

template <typename K, typename T, typename HF, typename EQF>
typename xmt_hash_table<K, T, HF, EQF>::size_type
xmt_hash_table<K, T, HF, EQF>::num_deleted_entries() const
{
  size_type num_deleted = 0;
  for (size_type i = 0; i < table_size; ++i)
  {
    num_deleted += (occupied[i] == XMT_HT_DELETED);
  }
  return num_deleted;
}

/***/

template <typename K, typename T, typename HF, typename EQF>
T&
xmt_hash_table<K, T, HF, EQF>::operator[](const K& key)
{
  // Find the element.
  size_type index = hash(key);
  size_type i = index;

  do
  {
    if (occupied[i] == XMT_HT_OCCUPIED)
    {
      if (compare(table[i].first, key)) return table[i].second;
    }
    else if (occupied[i] == XMT_HT_EMPTY)
    {
      // Acquire the lock.
      long probed = mt_readfe(occupied[i]);

      // Make sure the entry is still empty.
      if (probed == XMT_HT_EMPTY)
      {
        // Add the element at the current position in the table.
        const_cast<K&>(table[i].first) = key;
        table[i].second = T();
        mt_write(occupied[i], XMT_HT_OCCUPIED);

        return table[i].second;
      }
      else if (probed == XMT_HT_OCCUPIED)
      {
        if (compare(table[i].first, key))
        {
          mt_write(occupied[i], XMT_HT_OCCUPIED);
          return table[i].second;
        }
      }

      // The entry had been grabbed by someone else.  Release the lock.
      mt_write(occupied[i], probed);
    }

    i = (i + 1) & mask;
  } while (i != index);

  // TODO:  There is an error here.  If the element doesn't exist in the table
  //        and the table is full, a reference can't be returned.  In this
  //        case, we just need to die for now.  This would be fixed if
  //        dynamic reallocation were added.
  assert(false);

  // This line is only here to get rid of a compiler warning.  The code should
  // never get here.
  return table[i].second;
}

/***/

template <typename K, typename T, typename HF, typename EQF>
bool xmt_hash_table<K, T, HF, EQF>::member(const K& key) const
{
  size_type index = hash(key);
  size_type i = index;

  do
  {
    if (occupied[i] == XMT_HT_OCCUPIED)
    {
      if (compare(table[i].first, key)) return true;
    }
    else if (occupied[i] == XMT_HT_EMPTY)
    {
      break;
    }

    i = (i + 1) & mask;
  } while (i != index);

  return false;
}

/***/

template <typename K, typename T, typename HF, typename EQF>
bool xmt_hash_table<K, T, HF, EQF>::lookup(const K& key, T& value) const
{
  size_type index = hash(key);
  size_type i = index;

  do
  {
    if (occupied[i] == XMT_HT_OCCUPIED)
    {
      if (compare(table[i].first, key))
      {
        value = table[i].second;
        return true;
      }
    }
    else if (occupied[i] == XMT_HT_EMPTY)
    {
      break;
    }

    i = (i + 1) & mask;
  } while (i != index);

  return false;
}

/***/

template <typename K, typename T, typename HF, typename EQF>
pair<typename xmt_hash_table<K, T, HF, EQF>::value_type*, bool>
xmt_hash_table<K, T, HF, EQF>::insert(const K& key, const T& value)
{
  size_type index = hash(key);
  size_type i = index;

  do
  {
    // Entry is empty.  Add the element.
    if (occupied[i] == XMT_HT_EMPTY)
    {
      // Acquire the lock.
      long probed = mt_readfe(occupied[i]);

      // Make sure the entry is still empty.
      if (probed == XMT_HT_EMPTY)
      {
        // Add the element at the current position in the table.
        const_cast<K&>(table[i].first) = key;
        table[i].second = value;
        mt_write(occupied[i], XMT_HT_OCCUPIED);

        return pair<value_type*, bool>(table + i, true);
      }
      else
      {
        // If this entry is occupied and we match the key, exit without
        // inserting.
        if (probed == XMT_HT_OCCUPIED && compare(table[i].first, key))
        {
          mt_write(occupied[i], probed);
          return pair<value_type*, bool>(table + i, false);
        }
      }

      // The entry had been grabbed by someone else.  Release the lock.
      mt_write(occupied[i], probed);
    }
    else
    {
      // If this entry is occupied and we match the key, exit without
      // inserting.
      if (occupied[i] == XMT_HT_OCCUPIED && compare(table[i].first, key))
      {
        return pair<value_type*, bool>(table + i, false);
      }
    }

    i = (i + 1) & mask;
  } while (i != index);

  return pair<value_type*, bool>(NULL, false);
}

/***/

template <typename K, typename T, typename HF, typename EQF>
template <typename Vis>
pair<typename xmt_hash_table<K, T, HF, EQF>::value_type*, bool>
xmt_hash_table<K, T, HF, EQF>::insert(const K& key, const T& value, Vis visitor)
{
  size_type index = hash(key);
  size_type i = index;

  do
  {
    // Entry is empty.  Add the element.
    if (occupied[i] == XMT_HT_EMPTY)
    {
      // Acquire the lock.
      long probed = mt_readfe(occupied[i]);

      // Make sure the entry is still empty.
      if (probed == XMT_HT_EMPTY)
      {
        // Add the element at the current position in the table.
        const_cast<K&>(table[i].first) = key;
        table[i].second = value;
        visitor(table[i].second);
        mt_write(occupied[i], XMT_HT_OCCUPIED);

        return pair<value_type*, bool>(table + i, true);
      }
      else
      {
        // If this entry is occupied and we match the key, exit without
        // inserting.
        if (probed == XMT_HT_OCCUPIED && compare(table[i].first, key))
        {
          mt_write(occupied[i], probed);
          return pair<value_type*, bool>(table + i, false);
        }
      }

      // The entry had been grabbed by someone else.  Release the lock.
      mt_write(occupied[i], probed);
    }
    else
    {
      // If this entry is occupied and we match the key, exit without
      // inserting.
      if (occupied[i] == XMT_HT_OCCUPIED && compare(table[i].first, key))
      {
        return pair<value_type*, bool>(table + i, false);
      }
    }

    i = (i + 1) & mask;
  } while (i != index);

  return pair<value_type*, bool>(NULL, false);
}

/***/

template <typename K, typename T, typename HF, typename EQF>
pair<typename xmt_hash_table<K, T, HF, EQF>::value_type*, bool>
xmt_hash_table<K, T, HF, EQF>::insert_reuse(const K& key, const T& value)
{
  size_type index = hash(key);
  size_type i = index;
  size_type insert_pos = (std::numeric_limits<size_type>::max)();

  do
  {
    // Entry is empty.  Add the element.
    if (occupied[i] == XMT_HT_EMPTY)
    {
      // We found a deleted entry earlier.  Add the element at that position
      // in the table.
      if (insert_pos != (std::numeric_limits<size_type>::max)())
      {
        // Acquire the lock.
        long probed = mt_readfe(occupied[insert_pos]);

        // Make sure the entry is still deleted.
        if (probed == XMT_HT_DELETED)
        {
          // Add the element at the current position in the table.
          const_cast<K&>(table[insert_pos].first) = key;
          table[insert_pos].second = value;
          mt_write(occupied[insert_pos], XMT_HT_OCCUPIED);

          return pair<value_type*, bool>(table + insert_pos, true);
        }
        else if (compare(table[insert_pos].first, key))
        {
          // This entry is occupied, and we matched the key.  Exit without
          // inserting.
          mt_write(occupied[insert_pos], probed);
          return pair<value_type*, bool>(table + insert_pos, false);
        }

        // We need to start searching from insert_pos again.
        i = insert_pos;

        // The entry had been grabbed by someone else.  Release the lock.
        mt_write(occupied[insert_pos], probed);

        // Reset the insert pos back to max
        insert_pos = (std::numeric_limits<size_type>::max)();
      }
      else
      {
        // Acquire the lock.
        long probed = mt_readfe(occupied[i]);

        // Make sure the entry is still empty.
        if (probed == XMT_HT_EMPTY)
        {
          // Add the element at the current position in the table.
          const_cast<K&>(table[i].first) = key;
          table[i].second = value;
          mt_write(occupied[i], XMT_HT_OCCUPIED);

          return pair<value_type*, bool>(table + i, true);
        }
        else if (probed == XMT_HT_OCCUPIED)
        {
          // This entry is occupied.  If we match the key, exit without
          // inserting.
          if (compare(table[i].first, key))
          {
            mt_write(occupied[i], probed);
            return pair<value_type*, bool>(table + i, false);
          }
        }
        // We know that at this point the entry contains a deleted value.
        else if (insert_pos == (std::numeric_limits<size_type>::max)())
        {
          // This is the first deleted entry that has been encountered.  Keep
          // the lock here (because this is where the new element should be
          // inserted) and continue searching until either the element is found
          // or an empty entry is discovered.
          insert_pos = i;
        }

        // The entry had been grabbed by someone else.  Release the lock.
        mt_write(occupied[i], probed);
      }
    }
    else if (occupied[i] == XMT_HT_OCCUPIED)
    {
      // This entry is occupied.  If we match the key, exit without inserting.
      if (compare(table[i].first, key))
      {
        return pair<value_type*, bool>(table + i, false);
      }
    }
    // We know that at this point the entry contains a deleted value.
    else if (insert_pos == (std::numeric_limits<size_type>::max)())
    {
      // This is the first deleted entry that has been encountered.  Keep
      // the lock here (because this is where the new element should be
      // inserted) and continue searching until either the element is found
      // or an empty entry is discovered.
      insert_pos = i;
    }

    i = (i + 1) & mask;
  } while (i != index);

  // If there is still an insert_pos marked, then the table contained only
  // occupied and deleted entries.  Try to insert at the remembered deleted
  // location.
  if (insert_pos != (std::numeric_limits<size_type>::max)())
  {
    // Acquire the lock.
    long probed = mt_readfe(occupied[insert_pos]);

    // Make sure the entry is still deleted.
    if (probed == XMT_HT_DELETED)
    {
      // Add the element at the current position in the table.
      const_cast<K&>(table[insert_pos].first) = key;
      table[insert_pos].second = value;
      mt_write(occupied[insert_pos], XMT_HT_OCCUPIED);

      return pair<value_type*, bool>(table + insert_pos, true);
    }
  }

  return pair<value_type*, bool>(NULL, false);
}

/***/

template <typename K, typename T, typename HF, typename EQF>
template <typename Vis>
pair<typename xmt_hash_table<K, T, HF, EQF>::value_type*, bool>
xmt_hash_table<K, T, HF, EQF>::insert_reuse(const K& key, const T& value,
                                            Vis visitor)
{
  size_type index = hash(key);
  size_type i = index;
  size_type insert_pos = (std::numeric_limits<size_type>::max)();

  do
  {
    // Entry is empty.  Add the element.
    if (occupied[i] == XMT_HT_EMPTY)
    {
      // We found a deleted entry earlier.  Add the element at that position
      // in the table.
      if (insert_pos != (std::numeric_limits<size_type>::max)())
      {
        // Acquire the lock.
        long probed = mt_readfe(occupied[insert_pos]);

        // Make sure the entry is still deleted.
        if (probed == XMT_HT_DELETED)
        {
          // Add the element at the current position in the table.
          const_cast<K&>(table[insert_pos].first) = key;
          table[insert_pos].second = value;
          visitor(table[insert_pos].second);
          mt_write(occupied[insert_pos], XMT_HT_OCCUPIED);

          return pair<value_type*, bool>(table + insert_pos, true);
        }
        else if (compare(table[insert_pos].first, key))
        {
          // This entry is occupied, and we matched the key.  Exit without
          // inserting.
          mt_write(occupied[insert_pos], probed);
          return pair<value_type*, bool>(table + insert_pos, false);
        }

        // We need to start searching from insert_pos again.
        i = insert_pos;

        // The entry had been grabbed by someone else.  Release the lock.
        mt_write(occupied[insert_pos], probed);

        // Reset the insert pos back to max
        insert_pos = (std::numeric_limits<size_type>::max)();
      }
      else
      {
        // Acquire the lock.
        long probed = mt_readfe(occupied[i]);

        // Make sure the entry is still empty.
        if (probed == XMT_HT_EMPTY)
        {
          // Add the element at the current position in the table.
          const_cast<K&>(table[i].first) = key;
          table[i].second = value;
          visitor(table[i].second);
          mt_write(occupied[i], XMT_HT_OCCUPIED);

          return pair<value_type*, bool>(table + i, true);
        }
        else if (probed == XMT_HT_OCCUPIED)
        {
          // This entry is occupied.  If we match the key, exit without
          // inserting.
          if (compare(table[i].first, key))
          {
            mt_write(occupied[i], probed);
            return pair<value_type*, bool>(table + i, false);
          }
        }
        // We know that at this point the entry contains a deleted value.
        else if (insert_pos == (std::numeric_limits<size_type>::max)())
        {
          // This is the first deleted entry that has been encountered.  Keep
          // the lock here (because this is where the new element should be
          // inserted) and continue searching until either the element is found
          // or an empty entry is discovered.
          insert_pos = i;
        }

        // The entry had been grabbed by someone else.  Release the lock.
        mt_write(occupied[i], probed);
      }
    }
    else if (occupied[i] == XMT_HT_OCCUPIED)
    {
      // This entry is occupied.  If we match the key, exit without inserting.
      if (compare(table[i].first, key))
      {
        return pair<value_type*, bool>(table + i, false);
      }
    }
    // We know that at this point the entry contains a deleted value.
    else if (insert_pos == (std::numeric_limits<size_type>::max)())
    {
      // This is the first deleted entry that has been encountered.  Keep
      // the lock here (because this is where the new element should be
      // inserted) and continue searching until either the element is found
      // or an empty entry is discovered.
      insert_pos = i;
    }

    i = (i + 1) & mask;
  } while (i != index);

  // If there is still an insert_pos marked, then the table contained only
  // occupied and deleted entries.  Try to insert at the remembered deleted
  // location.
  if (insert_pos != (std::numeric_limits<size_type>::max)())
  {
    // Acquire the lock.
    long probed = mt_readfe(occupied[insert_pos]);

    // Make sure the entry is still deleted.
    if (probed == XMT_HT_DELETED)
    {
      // Add the element at the current position in the table.
      const_cast<K&>(table[insert_pos].first) = key;
      table[insert_pos].second = value;
      visitor(table[insert_pos].second);
      mt_write(occupied[insert_pos], XMT_HT_OCCUPIED);

      return pair<value_type*, bool>(table + insert_pos, true);
    }
  }

  return pair<value_type*, bool>(NULL, false);
}

/***/

template <typename K, typename T, typename HF, typename EQF>
bool xmt_hash_table<K, T, HF, EQF>::erase(const K& key)
{
  size_type index = hash(key);
  size_type i = index;

  do
  {
    // Entry is occupied.
    if (occupied[i] == XMT_HT_OCCUPIED)
    {
      // Match the entry.
      if (compare(table[i].first, key))
      {
        // Acquire the lock.
        long probed = mt_readfe(occupied[i]);

        // If the entry was deleted already, we don't care.  If the entry was
        // deleted and then the same key inserted again, we don't care.  Just
        // delete it again.  However, if it was deleted and then another key
        // inserted in the same entry, we need to not delete the new key.
        if (compare(table[i].first, key))
        {
          // Set the entry to deleted and release the lock.
          mt_write(occupied[i], XMT_HT_DELETED);

          return true;
        }

        // The entry had been deleted by another thread and then a different
        // key inserted by a second thread.  Release the lock.
        mt_write(occupied[i], probed);
      }
    }
    else if (occupied[i] == XMT_HT_EMPTY)
    {
      return false;
    }

    i = (i + 1) & mask;
  } while (i != index);

  return false;
}

/***/

template <typename K, typename T, typename HF, typename EQF>
bool xmt_hash_table<K, T, HF, EQF>::update(const K& key, const T& value)
{
  size_type index = hash(key);
  size_type i = index;

  do
  {
    if (occupied[i] == XMT_HT_OCCUPIED)
    {
      if (compare(table[i].first, key))
      {
        table[i].second = value;
        return true;
      }
    }
    else if (occupied[i] == XMT_HT_EMPTY)
    {
      break;
    }

    i = (i + 1) & mask;
  } while (i != index);

  return false;
}

/***/

template <typename K, typename T, typename HF, typename EQF>
template <typename Vis>
bool xmt_hash_table<K, T, HF, EQF>::update(const K& key, const T& value,
                                           Vis visitor)
{
  size_type index = hash(key);
  size_type i = index;

  do
  {
    if (occupied[i] == XMT_HT_OCCUPIED)
    {
      if (compare(table[i].first, key))
      {
        visitor(table[i].second, value);
        return true;
      }
    }
    else if (occupied[i] == XMT_HT_EMPTY)
    {
      break;
    }

    i = (i + 1) & mask;
  } while (i != index);

  return false;
}

/***/

template <typename K, typename T, typename HF, typename EQF>
bool xmt_hash_table<K, T, HF, EQF>::update_insert(const K& key, const T& value)
{
  size_type index = hash(key);
  size_type i = index;

  do
  {
    if (occupied[i] == XMT_HT_OCCUPIED)
    {
      if (compare(table[i].first, key))
      {
        table[i].second = value;
        return false;
      }
    }
    else if (occupied[i] == XMT_HT_EMPTY)
    {
      // Acquire the lock.
      long probed = mt_readfe(occupied[i]);

      // Make sure the entry is still empty.
      if (probed == XMT_HT_EMPTY)
      {
        // Add the element at the current position in the table.
        const_cast<K&>(table[i].first) = key;
        table[i].second = value;
        mt_write(occupied[i], XMT_HT_OCCUPIED);

        return true;
      }
      else if (probed == XMT_HT_OCCUPIED)
      {
        if (compare(table[i].first, key))
        {
          mt_write(occupied[i], XMT_HT_OCCUPIED);
          table[i].second = value;
          return false;
        }
      }

      // The entry had been grabbed by someone else.  Release the lock.
      mt_write(occupied[i], probed);
    }

    i = (i + 1) & mask;
  } while (i != index);

  return false;
}

/***/

template <typename K, typename T, typename HF, typename EQF>
template <typename Vis>
bool xmt_hash_table<K, T, HF, EQF>::update_insert(const K& key, const T& value,
                                                  Vis visitor)
{
  size_type index = hash(key);
  size_type i = index;

  do
  {
    if (occupied[i] == XMT_HT_OCCUPIED)
    {
      if (compare(table[i].first, key))
      {
        visitor(table[i].second, value);
        return false;
      }
    }
    else if (occupied[i] == XMT_HT_EMPTY)
    {
      // Acquire the lock.
      long probed = mt_readfe(occupied[i]);

      // Make sure the entry is still empty.
      if (probed == XMT_HT_EMPTY)
      {
        // Add the element at the current position in the table.
        const_cast<K&>(table[i].first) = key;
        table[i].second = value;
        mt_write(occupied[i], XMT_HT_OCCUPIED);

        return true;
      }
      else if (probed == XMT_HT_OCCUPIED)
      {
        if (compare(table[i].first, key))
        {
          mt_write(occupied[i], XMT_HT_OCCUPIED);
          visitor(table[i].second, value);
          return false;
        }
      }

      // The entry had been grabbed by someone else.  Release the lock.
      mt_write(occupied[i], probed);
    }

    i = (i + 1) & mask;
  } while (i != index);

  return false;
}

/***/

template <typename K, typename T, typename HF, typename EQF>
void
xmt_hash_table<K, T, HF, EQF>::swap(xmt_hash_table& rhs)
{
  size_type temp_mask = rhs.mask;
  value_type* temp_table = rhs.table;
  long* temp_occupied = rhs.occupied;
  size_type temp_table_size = rhs.table_size;

  rhs.mask = mask;
  rhs.table = table;
  rhs.occupied = occupied;
  rhs.table_size = table_size;

  mask = temp_mask;
  table = temp_table;
  occupied = temp_occupied;
  table_size = temp_table_size;
}

/***/

template <typename K, typename T, typename HF, typename EQF>
void xmt_hash_table<K, T, HF, EQF>::clear()
{
  memset(occupied, XMT_HT_EMPTY, table_size * sizeof(long));
}

/***/

template <typename K, typename T, typename HF, typename EQF>
template <typename Vis>
void xmt_hash_table<K, T, HF, EQF>::visit(Vis visitor)
{
  #pragma mta assert parallel
  for (size_type i = 0; i < table_size; ++i)
  {
    if (occupied[i] == XMT_HT_OCCUPIED) visitor(table + i);
  }
}

/***/

template <typename K, typename T, typename HF, typename EQF>
void xmt_hash_table<K, T, HF, EQF>::print() const
{
  std::cout << "Elements: " << size() << std::endl;
  for (size_type i = 0; i < table_size; ++i)
  {
    if (occupied[i] == XMT_HT_OCCUPIED)
    {
      std::cout << i << ": (" << table[i].first << ", " << table[i].second
                << ")" << std::endl;
    }
    else if (occupied[i] == XMT_HT_DELETED)
    {
      std::cout << i << ": Tombstoned." << std::endl;
    }
  }
}

/***/

template <typename K, typename T, typename HF, typename EQF>
unsigned long
xmt_hash_table<K, T, HF, EQF>::get_mmap_size()
{
  return 4 * sizeof(unsigned long) +
         table_size * (sizeof(value_type) + sizeof(long));
}

/***/

template <typename K, typename T, typename HF, typename EQF>
void xmt_hash_table<K, T, HF, EQF>::write_mmap(void* mapped_mem)
{
  // Can't mmap something that doesn't have a defined type.
  assert(mmap_traits<xmt_hash_table>::type);

  unsigned long* ul_mapped_mem = reinterpret_cast<unsigned long*>(mapped_mem);

  // Write the mmap type, mmap size, table_size, and mask to the mapped
  // memory.
  ul_mapped_mem[0] = mmap_traits<xmt_hash_table>::type;
  ul_mapped_mem[1] = get_mmap_size();
  ul_mapped_mem[2] = table_size;
  ul_mapped_mem[3] = mask;

  // Get pointers to the locations in the mapped memory for the arrays.
  value_type* table_ptr = reinterpret_cast<value_type*>(ul_mapped_mem + 4);
  long* occupied_ptr = reinterpret_cast<long*>(table_ptr + table_size);

  // Write the values for the arrays to the mapped memory.
  memcpy(table_ptr, table, table_size * sizeof(value_type));
  memcpy(occupied_ptr, occupied, table_size * sizeof(long));
}

/***/

template <typename K, typename T, typename HF, typename EQF>
void xmt_hash_table<K, T, HF, EQF>::read_mmap(void* mapped_mem)
{
  clear_mem();

  unsigned long* ul_mapped_mem = reinterpret_cast<unsigned long*>(mapped_mem);

  // Set the table size and mask.
  table_size = ul_mapped_mem[2];
  mask = ul_mapped_mem[3];

  // Get pointers to the locations in the mapped memory for the arrays.
  value_type* table_ptr = reinterpret_cast<value_type*>(ul_mapped_mem + 4);
  long* occupied_ptr = reinterpret_cast<long*>(table_ptr + table_size);

  // Set the pointers to the arrays to point to the mapped memory.
  table = table_ptr;
  occupied = occupied_ptr;
  initialized_by_mmap = true;
}

/***/

template <typename K, typename T, typename HF, typename EQF>
void
xmt_hash_table<K, T, HF, EQF>::deep_copy(const xmt_hash_table& rhs)
{
  #pragma mta assert par_newdelete
  table = new value_type[table_size];
  occupied = (long*) malloc(table_size * sizeof(long));

  size_type ts = table_size;

  #pragma mta assert nodep
  for (size_type i = 0; i < ts; ++i)
  {
    occupied[i] = rhs.occupied[i];
  }

  #pragma mta assert nodep
  for (size_type i = 0; i < ts; ++i)
  {
    if (occupied[i])
    {
      const_cast<K&>(table[i].first) = rhs.table[i].first;
      table[i].second = rhs.table[i].second;
    }
  }
}

/***/

template <typename K, typename T, typename HF, typename EQF>
void
xmt_hash_table<K, T, HF, EQF>::clear_mem()
{
  if (initialized_by_mmap)
  {
    initialized_by_mmap = false;
  }
  else
  {
    #pragma mta assert par_newdelete
    delete [] table;
    free(occupied);
  }
}

/***/

template <typename K, typename T, typename HF, typename EQF>
class mmap_traits<xmt_hash_table<K, T, HF, EQF> > {
public:
  // If either K or T is an undefined type, then the hash table's type should
  // also be undefined.
  static const unsigned long type =
    (mmap_traits<K>::type != MMAP_TYPE_NOT_DEFINED) *
    (mmap_traits<T>::type != MMAP_TYPE_NOT_DEFINED) *
    (20000 + mmap_traits<K>::type * 100 + mmap_traits<T>::type);
};

}

#undef XMT_HT_EMPTY
#undef XMT_HT_DELETED
#undef XMT_HT_OCCUPIED

#endif
