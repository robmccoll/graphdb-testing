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
/*! \file xmt_hash_set.hpp

    \brief A thread-safe set abstraction implemented using open address
           hashing.

    \author Greg Mackey (gemacke@sandia.gov)
    \author Eric Goodman (elgoodm@sandia.gov)
    \author Jon Berry (jberry@sandia.gov)

    \date 6/13/2008
*/
/****************************************************************************/

#ifndef MTGL_XMT_HASH_SET_HPP
#define MTGL_XMT_HASH_SET_HPP

#include <limits>
#include <iostream>

#include <mtgl/util.hpp>
#include <mtgl/partitioning.hpp>
#include <mtgl/hash_defs.hpp>

#define XMT_HASH_SET_INIT_SIZE 1024

namespace mtgl {

/***/

/*! \brief A set data structure implemented using hashing with open
           addressing.

    \author Greg Mackey (gemacke@sandia.gov)
    \author Eric Goodman (elgoodm@sandia.gov)
    \author Jon Berry (jberry@sandia.gov)

    \tparam K Key type.
    \tparam HF Function object type defining the hash function.
    \tparam EQF Function object type defining the equality function for keys.

    The hash set uses fine-grained locking at the entry level.  It also uses
    the double locking method developed by Eric Goodman of Sandia National
    Labs and David Haglin of Pacific Northwest National Labs.  Both locking
    methods increase the parallelism of the hash set.

    The set enforces unique keys, and it allows deletes.  Note that keys
    cannot be modified.  Once they are inserted into the hash, they can only
    be deleted.  Because performing an int_fetch_add during insert on an
    element counter causes a hotspot, the set doesn't keep track of the number
    of elements currently in itself.  Instead, a linear cost size() function
    is provided.

    The hash set performs deletions by marking the entries as deleted.  The
    normal insert function ignore these deleted entries when searching for an
    open spot in the table.  This method of deletion is known as tombstoning.
    Once a position in the table is deleted, the position is unavailable for
    use by later insertions.  Thus, erasing elements shrinks the size of the
    table.  This is done because reusing the elements makes the inserts
    more expensive in instruction count, loads, and stores, although the time
    performance was the same in my tests for inserts into an empty table.  If
    space is more of a concern than performance, an insert function,
    insert_reuse(), that reuses the deleted positions is provided.

    There are two main classes of functions: those designed to be called in a
    parallel context and those designed to be called in a serial context.  The
    parallel context functions can be further divided into functions that
    modify the keys and functions that don't.  Here is the division of
    functions:

      - Parallel Context Functions
         - Modify the Key
            - insert()
            - insert_reuse()
            - erase()

         - Don't Modify the Key
            - member()

      - Serial Context Functions
         - Constructors
         - Destructor
         - operator=()
         - size()
         - resize()
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
    function.  Most functions seem to get the best performance when block
    scheduled.  Here is an example of parallelizing insert, a parallel context
    function that modifies the key.

      #pragma mta block schedule
      #pragma mta assert parallel
      for (int i = 0; i < 1000; ++i)
      {
        xhs.insert(i);
      }

    Here is an example of parallelizing member, a parallel context function
    that doesn't modify the key.

      #pragma mta block schedule
      #pragma mta assert nodep
      for (int i = 0; i < 1000; ++i)
      {
        bool isMember = xhs.member(i);
      }

    The hash set provides two types of iteration: loop-based and stream-based.
    Loop-based iteration is provided by the visit function.  Visitors are
    used for loop-based parallelism because iterators for single element
    access cause considerable complication with thread safety.  The visit
    function loops over the table and applies a visitor function object to
    each occupied element in the table.  The visitor has access to a single
    element at a time.  The () operator of the visitor should have a single
    parameter that is a pointer to a value_type, which is the key.  Below is
    an example of a visitor object that finds the sum of the keys.

      template <typename HS>
      class set_visitor {
      public:
        set_visitor() : sum(0) {}

        void operator()(typename HS::value_type* k) { mt_incr(sum, *k); }

        int sum;
      };

    Now an example of using the visit function.

      xmt_hash_set<int> hs;
      <Items inserted into the hash here.>
      set_visitor<xmt_hash_set<int> > sv;
      hs.visit(sv);

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
    that takes a hash set and finds the sum of the keys.  This is a better way
    to find the sum of the keys than using the visit funciton because it
    avoids the hotspot of the increment on sum.

      template<typename HS>
      inline
      typename HS::value_type get_key_sum(HS& xhs)
      {
        typename HS::size_type stream_id = 0;
        typename HS::size_type num_streams = 1;

        typename HS::value_type sum = 0;

        #pragma mta for all streams stream_id of num_streams
        {
          typename HS::thread_iterator begin_iter =
              xhs.begin(stream_id, num_streams);

          typename HS::thread_iterator end_iter =
              xhs.end(stream_id, num_streams);

          typename HS::value_type my_sum = 0;

          // Find the sum of the keys for the occupied entries in this
          // stream's region.
          for ( ; begin_iter != end_iter; ++begin_iter) my_sum += *begin_iter;

          mt_incr(sum, my_sum);
        }

        return sum;
      }

    The hash set offers static resizing via the resize function.  Note that
    this function is the only way to resize the table.  Also, it is not
    thread-safe in the sense that other operations on the hash set cannot be
    called at the same time as the resize function.  The user must ensure that
    the resize function is not interleaved with other hash set functions.
    (i.e. The user shouldn't put the resize function inside a loop.)  Making
    the resize function completely thread-safe would require implementing a
    reader-writer locking scheme which can significantly slow down the hash
    set performance.
*/
template <typename K, typename HF = default_hash_func<K>,
          typename EQF = default_eqfcn<K> >
class xmt_hash_set {
public:
  typedef K key_type;
  typedef const K value_type;
  typedef EQF key_compare;
  typedef hash_size_type size_type;
  class thread_iterator;

  /// \brief The normal constructor.
  /// \param size The initial allocation size.
  inline xmt_hash_set(size_type size = XMT_HASH_SET_INIT_SIZE);

  /// \brief The copy constructor, which performs a deep copy.
  /// \param h The xmt_hash_set to be copied.
  inline xmt_hash_set(const xmt_hash_set& h);

  /// The destructor.
  inline ~xmt_hash_set();

  /// \brief The assignment operator, which performs a deep copy.
  /// \param rhs The xmt_hash_set to be copied.
  inline xmt_hash_set& operator=(const xmt_hash_set& rhs);

  /// Returns the beginning iterator for this stream.
  inline thread_iterator begin(size_type stream_id, size_type num_streams);

  /// Returns the ending iterator for this stream.
  inline thread_iterator end(size_type stream_id, size_type num_streams);

  /// \brief Returns the number of elements in the hash set.
  /// \return Number of elements in the hash set.
  ///
  /// This is a linear cost operation.
  inline size_type size() const;

  /// \brief Resizes the hash set to have new_size slots.  If there are
  ///        any entries in the set, they are rehased to the new set.
  /// \param new_size The new size for the hash set.
  inline void resize(size_type new_size);

  /// \brief Returns the allocation size of the set.
  /// \return Allocation size of set.
  size_type capacity() const { return table_size; }

  /// \brief Returns true if there are no elements in the hash set.
  /// \return True if the set is empty.  Otherwise, false.
  ///
  /// This is a linear cost operation.
  inline size_type empty() const { return size() == 0; }

  /// \brief Tests if key exists in the set.
  /// \param key The search key.
  /// \return True if the key exists in the set.  Otherwise, false.
  inline bool member(const key_type& key) const;

  /// \brief Inserts an element into the hash set.
  /// \param key The insertion key.
  /// \return Pair of value_type* and bool.  The first is a pointer to the
  ///         inserted key.  If the key already existed, it points to the
  ///         existing key.  If the table was full, it is set to NULL.  The
  ///         second is true if the key was inserted and false if it already
  ///         existed.
  inline pair<value_type*, bool> insert(const key_type& key);

  /// \brief Inserts an element into the hash set reusing deleted entries.
  /// \param key The insertion key.
  /// \return Pair of value_type* and bool.  The first is a pointer to the
  ///         inserted key.  If the key already existed, it points to the
  ///         existing key.  If the table was full, it is set to NULL.  The
  ///         second is true if the key was inserted and false if it already
  ///         existed.
  ///
  /// This version of insert is slightly more expensive than insert().
  inline pair<value_type*, bool> insert_reuse(const key_type& key);

  /// \brief Erases key from the set.
  /// \param key The delete key.
  /// \return True if the key existed in the set and was erased.  Otherwise,
  ///         false.
  inline bool erase(const key_type& key);

  /// Swaps the contents of the two hash sets.
  inline void swap(xmt_hash_set& rhs);

  /// Clears all entries from the hash set.
  inline void clear();

  /// \brief Apply the method encapsulated by "visitor" to each set element
  ///        in parallel.
  /// \param visitor A closure object containing a function to be applied to
  ///        each element of the set.
  template <typename Vis> inline void visit(Vis& visitor) const;

  /// Prints the keys in the hash set.
  void print() const;

private:
  /// \brief Performs a deep copy from rhs to the calling object.
  /// \param rhs The xmt_hash_set to be copied.
  inline void deep_copy(const xmt_hash_set& rhs);

  /// \brief Returns the index into table where key belongs.
  ///
  /// Ands the user supplied hash function with the mask to get an index
  /// within the range 0 <= index < capacity().
  size_type hash(const key_type& key) const { return hash_func(key) & mask; }

private:
  size_type mask;
  key_type* table;
  long* occupied;
  size_type table_size;
  EQF compare;
  HF hash_func;
};

/// \brief This is a forward iterator that will always be given in pairs to
///        represent a range of values in the hash set.  This iterator
///        should only be used in "for all streams i of n" sections of code.
///
/// This iterator will always point to an occupied entry in the set.  On
/// initialization, the iterator is advanced to the next occupied entry.  If
/// the end of the set is reached, the iterator will point to one past the
/// end of the set.
template <typename K,typename HF, typename EQF>
class xmt_hash_set<K, HF, EQF>::thread_iterator {
public:
  thread_iterator() : occupied(0), table(0), table_size(0), index(0) {}

  thread_iterator(xmt_hash_set& s, size_type idx) :
    occupied(s.occupied), table(s.table),
    table_size(s.table_size), index(idx)
  {
    while (index != table_size && occupied[index] != 2) ++index;
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
    } while (index != table_size && occupied[index] != 2);

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
    } while (index != table_size && occupied[index] != 2);

    return temp;
  }

  value_type& operator*() const { return table[index]; }
  value_type* operator->() const { return &(table[index]); }

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
      if (occupied[i] == 2) ++num_entries;
    }

    return num_entries;
  }

private:
  long* occupied;
  key_type* table;
  size_type table_size;
  size_type index;
};

/***/

template <typename K, typename HF, typename EQF>
xmt_hash_set<K, HF, EQF>::xmt_hash_set(size_type size)
{
  size_type value = 16;
  if (size < 16)  size = 16;
  while (value < size)  { value *= 2; }
  table_size = value;
  mask = value - 1;

  // Allocate memory for table and occupied vector, and initialize occupied
  // vector to all 0's (empty).
  table = (key_type*) malloc(table_size * sizeof(key_type));
  occupied = (long*) malloc(table_size * sizeof(long));

  #pragma mta assert parallel
  #pragma mta block schedule
  for (size_type i = 0; i < table_size; ++i) occupied[i] = 0;
}

/***/

template <typename K, typename HF, typename EQF>
typename xmt_hash_set<K, HF, EQF>::thread_iterator
xmt_hash_set<K, HF, EQF>::begin(size_type stream_id, size_type num_streams)
{
  size_type index = begin_block_range(table_size, stream_id, num_streams);
  return thread_iterator(*this, index);
}

/***/

template <typename K, typename HF, typename EQF>
typename xmt_hash_set<K, HF, EQF>::thread_iterator
xmt_hash_set<K, HF, EQF>::end(size_type stream_id, size_type num_streams)
{
  size_type index = end_block_range(table_size, stream_id, num_streams);
  return thread_iterator(*this, index);
}

/***/

template <typename K, typename HF, typename EQF>
xmt_hash_set<K, HF, EQF>::xmt_hash_set(const xmt_hash_set& h) :
    mask(h.mask), table_size(h.table_size)
{
  deep_copy(h);
}

/***/

template <typename K, typename HF, typename EQF>
xmt_hash_set<K, HF, EQF>::~xmt_hash_set()
{
  free(table);
  free(occupied);
}

/***/

template <typename K, typename HF, typename EQF>
xmt_hash_set<K, HF, EQF>&
xmt_hash_set<K, HF, EQF>::operator=(const xmt_hash_set& rhs)
{
  free(table);
  free(occupied);

  mask = rhs.mask;
  table_size = rhs.table_size;

  deep_copy(rhs);

  return *this;
}

/***/

template <typename K, typename HF, typename EQF>
typename xmt_hash_set<K, HF, EQF>::size_type
xmt_hash_set<K, HF, EQF>::size() const
{
  size_type num_elements = 0;
  for (size_type i = 0; i < table_size; ++i) num_elements += (occupied[i] == 2);
  return num_elements;
}

/***/

template <typename K, typename HF, typename EQF>
void
xmt_hash_set<K, HF, EQF>::resize(size_type new_size)
{
  // Get copy of old values.
  key_type* old_table = table;
  long* old_occupied = occupied;
  size_type old_table_size = table_size;

  size_type num_elements = size();

  // Find the new size of the set.
  size_type value = 16;
  if (new_size < 16)  new_size = 16;
  while (value < new_size || value < num_elements)  { value *= 2; }
  table_size = value;
  mask = value - 1;

  // Allocate memory for table and occupied vector, and initialize occupied
  // vector to all 0's (empty).
  table = (key_type*) malloc(table_size * sizeof(key_type));
  occupied = (long*) malloc(table_size * sizeof(long));

  #pragma mta assert parallel
  #pragma mta block schedule
  for (size_type i = 0; i < table_size; ++i) occupied[i] = 0;

  // If the set isn't empty, copy the existing entries.
  if (num_elements > 0)
  {
    // Insert all the existing entries into the new set.
    #pragma mta assert parallel
    #pragma mta block schedule
    for (size_type i = 0; i < old_table_size; ++i)
    {
      if (old_occupied[i] == 2) insert(old_table[i]);
    }
  }

  // Delete memory for the old table and occupied vector.
  free(old_table);
  free(old_occupied);
}

/***/

template <typename K, typename HF, typename EQF>
bool xmt_hash_set<K, HF, EQF>::member(const K& key) const
{
  size_type index = hash(key);
  size_type i = index;

  do
  {
    if (occupied[i] == 2)
    {
      if (compare(table[i], key)) return true;
    }
    else if (occupied[i] == 0)
    {
      break;
    }

    i = (i + 1) & mask;
  } while (i != index);

  return false;
}

/***/

template <typename K, typename HF, typename EQF>
pair<typename xmt_hash_set<K, HF, EQF>::value_type*, bool>
xmt_hash_set<K, HF, EQF>::insert(const K& key)
{
  size_type index = hash(key);
  size_type i = index;

  do
  {
    // Entry is empty.  Add the element.
    if (occupied[i] == 0)
    {
      // Acquire the lock.
      long probed = mt_readfe(occupied[i]);

      // Make sure the entry is still empty.
      if (probed == 0)
      {
        // Add the element at the current position in table.
        const_cast<K&>(table[i]) = key;
        mt_write(occupied[i], 2);

        return pair<value_type*, bool>(table + i, true);
      }
      else
      {
        // If this entry is occupied and we match the key, exit without
        // inserting.
        if (probed == 2 && compare(table[i], key))
        {
          mt_write(occupied[i], probed);
          return pair<value_type*, bool>(table + i, false);
        }
      }

      // The entry had been grabbed by someone else.  Release the lock.
      mt_write(occupied[i], probed);
    }
    // This entry is occupied.
    else
    {
      // We matched the key.  Exit without inserting.
      if (occupied[i] == 2 && compare(table[i], key))
      {
        return pair<value_type*, bool>(table + i, false);
      }
    }

    i = (i + 1) & mask;
  } while (i != index);

  return pair<value_type*, bool>(NULL, false);
}

/***/

template <typename K, typename HF, typename EQF>
pair<typename xmt_hash_set<K, HF, EQF>::value_type*, bool>
xmt_hash_set<K, HF, EQF>::insert_reuse(const K& key)
{
  size_type index = hash(key);
  size_type i = index;
  size_type insert_pos = (std::numeric_limits<size_type>::max)();

  do
  {
    // Entry is empty.  Add the element.
    if (occupied[i] == 0)
    {
      // We found a deleted entry earlier.  Add the element at that position
      // in table.
      if (insert_pos != (std::numeric_limits<size_type>::max)())
      {
        // Acquire the lock.
        long probed = mt_readfe(occupied[insert_pos]);

        // Make sure the entry is still deleted.
        if (probed == 1)
        {
          // Add the element at the current position in table.
          const_cast<K&>(table[insert_pos] = key);
          mt_write(occupied[insert_pos], 2);

          return pair<value_type*, bool>(table + insert_pos, true);
        }
        else if (compare(table[insert_pos], key))
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
      }
      else
      {
        // Acquire the lock.
        long probed = mt_readfe(occupied[i]);

        // Make sure the entry is still empty.
        if (probed == 0)
        {
          // Add the element at the current position in table.
          const_cast<K&>(table[i] = key);
          mt_write(occupied[i], 2);

          return pair<value_type*, bool>(table + i, true);
        }
        else if (probed == 2)
        {
          // This entry is occupied.  If we match the key, exit without
          // inserting.
          if (compare(table[i], key))
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
        
        // Reset the insert pos back to max
        insert_pos = (std::numeric_limits<size_type>::max)();
      }
    }
    // This entry is occupied.
    else if (occupied[i] == 2)
    {
      // We matched the key.  Exit without inserting.
      if (compare(table[i], key))
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

  return pair<value_type*, bool>(NULL, false);
}

/***/

template <typename K, typename HF, typename EQF>
bool xmt_hash_set<K, HF, EQF>::erase(const K& key)
{
  size_type index = hash(key);
  size_type i = index;

  do
  {
    // Entry is full.
    if (occupied[i] == 2)
    {
      // Match the entry.
      if (compare(table[i], key))
      {
        // Acquire the lock.
        long probed = mt_readfe(occupied[i]);

        // If the entry was deleted already, we don't care.  If the entry was
        // deleted and then the same key inserted again, we don't care.  Just
        // delete it again.  However, if it was deleted and then another key
        // inserted in the same entry, we need to not delete the new key.
        if (compare(table[i], key))
        {
          // Set the entry to deleted and release the lock.
          mt_write(occupied[i], 1);

          return true;
        }

        // The entry had been deleted by another thread and then a different
        // key inserted by a second thread.  Release the lock.
        mt_write(occupied[i], probed);
      }
    }
    else if (occupied[i] == 0)
    {
      return false;
    }

    i = (i + 1) & mask;
  } while (i != index);

  return false;
}

/***/

template <typename K, typename HF, typename EQF>
void
xmt_hash_set<K, HF, EQF>::swap(xmt_hash_set& rhs)
{
  size_type temp_mask = rhs.mask;
  key_type* temp_table = rhs.table;
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

template <typename K, typename HF, typename EQF>
void xmt_hash_set<K, HF, EQF>::clear()
{
  memset(occupied, 0, table_size * sizeof(long));
}

/***/

template <typename K, typename HF, typename EQF>
template <typename Vis>
void xmt_hash_set<K, HF, EQF>::visit(Vis& visitor) const
{
  #pragma mta assert parallel
  for (size_type i = 0; i < table_size; ++i)
  {
    if (occupied[i] == 2) visitor(table + i);
  }
}

/***/

template <typename K, typename HF, typename EQF>
void xmt_hash_set<K, HF, EQF>::print() const
{
  std::cout << "Elements: " << size() << std::endl;
  for (size_type i = 0; i < table_size; ++i)
  {
    if (occupied[i] == 2)
    {
      std::cout << i << ": " << table[i] << std::endl;
    }
    else if (occupied[i] == 1)
    {
      std::cout << i << ": Tombstoned." << std::endl;
    }
  }
}

/***/

template <typename K, typename HF, typename EQF>
void
xmt_hash_set<K, HF, EQF>::deep_copy(const xmt_hash_set& rhs)
{
  table = (key_type*) malloc(table_size * sizeof(key_type));
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
    if (occupied[i]) const_cast<K&>(table[i]) = rhs.table[i];
  }
}

/// \brief Puts the keys from the hash table in a contiguous array, and
///        returns the number of keys.
///
/// \param ht The hash table.
/// \param keys The keys array that gets overwritten.  This is assumed to be
///             size() or bigger.
/// \return The number of keys.
template <typename K, typename HF, typename EQF>
inline
typename xmt_hash_set<K, HF, EQF>::size_type
get_keys(xmt_hash_set<K, HF, EQF>& hs,
         typename xmt_hash_set<K, HF, EQF>::key_type* keys)
{
  typedef xmt_hash_set<K, HF, EQF> HS;
  typename HS::size_type num_keys = 0;

  typename HS::size_type stream_id = 0;
  typename HS::size_type num_streams = 1;

  #pragma mta use 100 streams
  #pragma mta for all streams stream_id of num_streams
  {
    typename HS::thread_iterator begin = hs.begin(stream_id, num_streams);
    typename HS::thread_iterator end = hs.end(stream_id, num_streams);

    // Claim a block of indices for the occupied entries in the stream's
    // region.
    typename HS::size_type start_index = mt_incr(num_keys, end - begin);

    // Get the keys from this stream's region.
    for ( ; begin != end; ++begin)
    {
      keys[start_index] = *begin;
      ++start_index;
    }
  }

  return num_keys;
}

}

#undef XMT_HASH_SET_INIT_SIZE

#endif
