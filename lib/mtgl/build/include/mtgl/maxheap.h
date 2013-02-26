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

// //////////////////////////////////////////////////////////////////////
// --- COPYRIGHT NOTICE ---------------------------------------------
// FastCommunityMH - infers community structure of networks
// Copyright (C) 2004 Aaron Clauset
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// See http://www.gnu.org/licenses/gpl.txt for more details.
//
// //////////////////////////////////////////////////////////////////////
// Author       : Aaron Clauset  (aaron@cs.unm.edu)                    //
// Location     : U. Michigan, U. New Mexico                           //
// Time         : January-August 2004                                  //
// Collaborators: Dr. Cris Moore (moore@cs.unm.edu)                    //
//                Dr. Mark Newman (mejn@umich.edu)                     //
// //////////////////////////////////////////////////////////////////////

/****************************************************************************/
/*! \file maxheap.h

    \brief This class is only here because it is needed by wcnm().  It
           should not be used for any other purpose.
*/
/****************************************************************************/

namespace mtgl {

extern double* a;         // A_i (for computing WT consolidation ratios)
                          // a is define in wcnm.hpp.

extern double* num_nodes; // (For computing WT consolidation ratios.)
                          // The number of nodes in each community.

extern double W;          // Sum of edge weights.

#if !defined(TUPLE_INCLUDED)
#define TUPLE_INCLUDED
struct tuple_o {
  double m;              // Stored value  (dQ for CNM).
  double sv_diff;        // Stored support variance incr.
  int i;                 // Row index.
  int j;                 // Column index.
  int k;                 // Heap index.

  double WT_obj;         // Wakita/Tsurumi objective
                         //   m * min(a_i / a_j, a_j / a_i)
                         //   - or -
                         //   m * min(num_nodes[i] / num_nodes[j],
                         //           num_nodes[j] / num_nodes[i])

  double MB_obj;         // McCloskey/Bader objective:
                         //   R_ij / std(R_ij), where R_ij = W * dQ
};
#endif

#if !defined(MAXHEAP_INCLUDED)
#define MAXHEAP_INCLUDED

#include <iostream>

struct hnode { tuple_o* d; };
const int heapmin = 3;
// const double tiny = -4294967296.0;

/*!
   Because using this heap requires us to be able to modify an arbitrary
   element's data in constant O(1) time, I use to tricky tactic of having
   elements in an array-based heap only contain addresses to the data, rather
   than the data itself. In this manner, an external program may freely modify
   an element of the heap, given that it possesses a pointer to said element
   (in an array-based heap, the addresses and the value in that address are
   not bound and thus may change during the heapify() operation).
*/
class maxheap {
private:
  hnode* A;              ///< Maxheap array
  int heaplimit;         ///< First unused element of heap array
  int arraysize;         ///< Size of array
  bool isempty;          ///< T if heap is empty; F otherwise

  int downsift(int i);   ///< Sift A[i] down in heap.
  int upsift(int i);     ///< Sift A[i] up in heap.
  int left(int i);       ///< Returns index of left child.
  int right(int i);      ///< Returns index of right child.
  int parent(int i);     ///< Returns index of parent.
  void grow();           ///< Increase size of array A.
  void shrink();         ///< Decrease size of array A.

public:
  /// Default constructor.
  maxheap();

  /// Main constructor.
  maxheap(int size);

  /// Default destructor.
  ~maxheap();

  /// Returns heaplimit value.
  int heapSize();

  /// Returns isempty value.
  bool heapIsEmpty();

  /// Heap-inserts newData, returns the address of it.
  tuple_o* insertItem(const tuple_o newData);

  /// Removes and returns heap max, reheapifies.
  tuple_o popMaximum();

  /// Returns the heap max; no change to heap.
  tuple_o returnMaximum();

  /// Displays contents of the heap.
  void printHeap();

  /// Displays top 10 entries in the heap.
  void printHeapTop10();

  /// Updates the value of the tuple at address.
  void updateItem(tuple_o* address, tuple_o newData);

  /// Update only the stored value of tuple at address.
  void updateItem(tuple_o* address, double newStored);

  /// Remove an item from the heap.
  void deleteItem(tuple_o* address);

  int returnArraysize() { return arraysize; }
  int returnHeaplimit() { return heaplimit; }
};

maxheap::maxheap()
{
  heaplimit = 1;              // First free location is A[1].
  arraysize = heapmin + 1;
  isempty = true;             // Heap is initially empty.
  A = new hnode[arraysize];   // Allocate array for heap.

  // Initialize heap values.
  for (int i = 0; i < arraysize; i++)
  {
    tuple_o* newtuple = new tuple_o;
    A[i].d = newtuple;

    // A very negative value; unlikely to be a valid dQ.
    A[i].d->m = -4294967296.0;

    A[i].d->i = 0;
    A[i].d->j = 0;
    A[i].d->k = i;
    A[i].d->WT_obj = A[i].d->m;
    A[i].d->MB_obj = A[i].d->m;
  }
}

maxheap::maxheap(int size)
{
  heaplimit = 1;              // First free location is A[1].
  arraysize = size + 1;
  isempty = true;             // Heap is initially empty.
  A = new hnode[arraysize];   // Allocate array for heap.

  // Initialize heap values.
  for (int i = 0; i < arraysize; i++)
  {
    tuple_o* newtuple = new tuple_o;
    A[i].d = newtuple;

    // A very negative value; unlikely to be a valid dQ.
    A[i].d->m = -4294967296.0;

    A[i].d->i = 0;
    A[i].d->j = 0;
    A[i].d->k = i;
    A[i].d->WT_obj = A[i].d->m;
    A[i].d->MB_obj = A[i].d->m;
  }
}

maxheap::~maxheap()
{
  // First delete the containers pointed to.
  for (int i = 0; i < arraysize; i++) delete A[i].d;

  // Then delete the list of pointers to containers.
  delete [] A;
}

tuple_o maxheap::popMaximum()                           // O(log k) time
{
  tuple_o temp = returnMaximum();
  deleteItem(A[1].d);
  return temp;
}

tuple_o maxheap::returnMaximum()                        // O(1) time
{
  tuple_o temp;

  // Grab A's data.
  temp.m = A[1].d->m;
  temp.i = A[1].d->i;
  temp.j = A[1].d->j;
  temp.k = A[1].d->k;
  temp.WT_obj = A[1].d->WT_obj;
  temp.MB_obj = A[1].d->MB_obj;

  return temp;
}

int maxheap::downsift(int index)                        // O(log k) time
{
  bool stopFlag = false;
  int L = left(index);
  int R = right(index);
  int swap;

  while (!stopFlag)
  {
    // Check that both children are within the array boundaries.
    if ((L < heaplimit) && (R < heaplimit))
    {
      // First, choose larger of the children.
#if defined(WT)
      if (A[L].d->WT_obj > A[R].d->WT_obj)
#elif defined(MB)
      if (A[L].d->MB_obj > A[R].d->MB_obj)
#else
      if (A[L].d->m > A[R].d->m)
#endif
      {
        swap = L;
      }
      else
      {
        swap = R;
      }
    }
    else
    {
      // Only one child to consider.
      if (L < heaplimit)
      {
        swap = L;
      }
      else
      {
        break;
      }
    }

   // Now decide if need to exchange A[index] with A[swap].
#if defined(WT)
    if (A[index].d->WT_obj < A[swap].d->WT_obj)
#elif defined(MB)
    if (A[index].d->MB_obj < A[swap].d->MB_obj)
#else
    if (A[index].d->m < A[swap].d->m)
#endif
    {
      // Exchange pointers A[index] and A[swap].
      tuple_o* temp = A[index].d;
      A[index].d = A[swap].d;
      A[index].d->k = index;    // Note A[index].d's change of array location.
      A[swap].d = temp;
      A[swap].d->k = swap;      // Note A[swap].d's change in array location.

      // Update indices for next pass.
      index = swap;
      L = left(index);
      R = right(index);
    }
    else
    {
      stopFlag = true;
    }
  }

  // Return the new index location of downsifted element.
  return index;
}

int maxheap::upsift(int index)                          // O(log k) time
{
  bool stopFlag = false;
  int P = parent(index);

  while (!stopFlag)
  {
    // Decide if A[index] needs to move up in tree.
#if defined(WT)
    if ((P > 0) && (A[index].d->WT_obj > A[P].d->WT_obj))
#elif defined(MB)
    if ((P > 0) && (A[index].d->MB_obj > A[P].d->MB_obj))
#else
    if ((P > 0) && (A[index].d->m > A[P].d->m))
#endif
    {
      // Exchange A[index] and A[P].
      tuple_o* temp = A[index].d;
      A[index].d = A[P].d;
      A[index].d->k = index;    // Note A[index].d's change of array location.
      A[P].d = temp;
      A[P].d->k = P;            // Note A[P].d's change of array location.

      // Update indices for next pass.
      index = P;
      P = parent(index);
    }
    else
    {
      stopFlag = true;
    }
  }

  return index;
}

int maxheap::left(int index) { return 2 * index; }
int maxheap::right(int index) { return 2 * index + 1; }
int maxheap::parent(int index) { return (int) index / 2; }

void maxheap::grow()                                    // O(k) time
{
  // Scratch space for expansion of A.
  hnode* B = new hnode[arraysize];

  // Copy A into B.
  for (int i = 0; i < arraysize; i++) B[i].d = A[i].d;

  // Delete old array of addresses.
  delete [] A;

  // Grow A by factor of 2.
  A = new hnode[2 * arraysize];

  // Copy B into first half of A.
  for (int i = 0; i < arraysize; i++) A[i].d = B[i].d;

  // Initialize new heap values.
  for (int i = arraysize; i < (2 * arraysize); i++)
  {
    tuple_o* newtuple = new tuple_o;
    A[i].d = newtuple;
    A[i].d->m = -4294967296.0;
    A[i].d->i = 0;
    A[i].d->j = 0;
    A[i].d->k = i;
    A[i].d->WT_obj = A[i].d->m;
    A[i].d->MB_obj = A[i].d->m;
  }

  // Delete scratch space B.
  delete [] B;

  // Update size of array.
  arraysize = 2 * arraysize;
}

void maxheap::shrink()                                  // O(k) time
{
  // Scratch space for contraction of A.
  hnode* B = new hnode[heaplimit];

  // Copy A into B.
  for (int i = 0; i < heaplimit; i++) B[i].d = A[i].d;

  // Delete the unused containers in A.
  for (int i = heaplimit; i < arraysize; i++) delete A[i].d;

  // Delete old array of addresses.
  delete [] A;

  // Shrink A by factor of 2.
  A = new hnode [arraysize / 2];

  // Copy B into A.
  for (int i = 0; i < heaplimit; i++) A[i].d = B[i].d;

  // Initialize new heap values.
  for (int i = heaplimit; i < (arraysize / 2); i++)
  {
    tuple_o* newtuple = new tuple_o;
    A[i].d = newtuple;
    A[i].d->m = -4294967296.0;
    A[i].d->i = 0;
    A[i].d->j = 0;
    A[i].d->k = i;
    A[i].d->WT_obj = A[i].d->m;
    A[i].d->MB_obj = A[i].d->m;
  }

  // Delete scratch space B.
  delete [] B;

  // Update size of array.
  arraysize = arraysize / 2;
}

tuple_o* maxheap::insertItem(const tuple_o newData)     // O(log k) time
{
  // If heap is full, grow by factor of 2.
  if (heaplimit >= (arraysize - 1)) grow();

  int index = heaplimit;

  // Copy newData onto the bottom of the heap.
  A[index].d->m = newData.m;
  A[index].d->i = newData.i;
  A[index].d->j = newData.j;
  A[index].d->k = index;

  double i_nodes = num_nodes[newData.i];
  double j_nodes = num_nodes[newData.j];
  double r1 = (i_nodes / j_nodes);
  double r2 = (j_nodes / i_nodes);

#ifdef MB
  double ai = a[newData.i];
  double aj = a[newData.j];
  double mij = 2 * W * ai * aj;
  double xij = ((newData.m + 2 * ai * aj)) * W;
  double std_rij = sqrt(mij);

  A[index].d->MB_obj = (xij - mij) / std_rij;
#endif

  A[index].d->WT_obj = newData.m * (std::min)(r1 * r1 * r1 * r1 * r1,
                                              r2 * r2 * r2 * r2 * r2);

  // Store pointer to container.
  tuple_o* pointer = A[index].d;

  // Note the larger heap.
  heaplimit++;

  // Upsift new item up to proper spot.
  upsift(index);

  if (heaplimit > 1) isempty = false;

  return pointer;
}

void maxheap::updateItem(tuple_o* address, tuple_o newData)  // O(log k) time
{
#if defined(WT)
  double oldval = address->WT_obj;
#elif defined(MB)
  double oldval = address->MB_obj;
#else
  double oldval = address->m;
#endif

  // Udpate the dQ stored.
  address->m = newData.m;

  // Udpate the community to which to join.
  address->j = newData.j;

#if defined(WT)
  double i_nodes = num_nodes[newData.i];
  double j_nodes = num_nodes[newData.j];
  double r1 = (i_nodes / j_nodes);
  double r2 = (j_nodes / i_nodes);

  address->WT_obj = newData.m * (std::min)(r1 * r1 * r1 * r1 * r1,
                                           r2 * r2 * r2 * r2 * r2);

  double newval = address->WT_obj;
#elif defined(MB)
  double ai = a[newData.i];
  double aj = a[newData.j];
  double mij = 2 * W * ai * aj;
  double xij = ((address->m + 2 * ai * aj)) * W;
  double std_rij = sqrt(mij);

  address->MB_obj = (xij - mij) / std_rij;

  double newval = address->MB_obj;
#else
  double newval = address->m;
#endif

  if (oldval > newval)
  {
    // Downsift if old value was larger.
    downsift(address->k);
  }
  else
  {
    // Upsift otherwise.
    upsift(address->k);
  }
}

void maxheap::updateItem(tuple_o* address, double newStored)  // O(log k) time
{
#if defined(WT)
  double oldval = address->WT_obj;
#elif defined(MB)
  double oldval = address->MB_obj;
#else
  double oldval = address->m;
#endif

  // Udpate the dQ stored.
  address->m = newStored;

#if defined(WT)
  double i_nodes = num_nodes[address->i];
  double j_nodes = num_nodes[address->j];
  double r1 = (i_nodes / j_nodes);
  double r2 = (j_nodes / i_nodes);

  address->WT_obj = newStored * (std::min)(r1 * r1 * r1 * r1 * r1,
                                           r2 * r2 * r2 * r2 * r2);

  double newval = address->WT_obj;
#elif defined(MB)
  double ai = a[address->i];
  double aj = a[address->j];
  double mij = 2 * W * ai * aj;
  double xij = ((newStored + 2 * ai * aj)) * W;
  double std_rij = sqrt(mij);

  address->MB_obj = (xij - mij) / std_rij;

  double newval = address->MB_obj;
#else
  double newval = address->m;
#endif

  if (oldval > newval)
  {
    // Downsift if old value was larger.
    downsift(address->k);
  }
  else
  {
    // Upsift otherwise.
    upsift(address->k);
  }
}

void maxheap::deleteItem(tuple_o* address)
{
  int small_val = heaplimit - 1;
  int index = address->k;

  // Check if deleting last item in heap.
  if (heaplimit == 2)
  {
    // Zero out root data.
    A[1].d->m = -4294967296.0;
    A[1].d->i = 0;
    A[1].d->j = 0;
    A[1].d->k = 1;
    A[1].d->WT_obj = -4294967296.0;
    A[1].d->MB_obj = -4294967296.0;

    // Note the heap's emptiness.
    isempty = true;

    // Decrement size of heap to empty.
    heaplimit--;
  }
  else
  {
    // Zero out the deleted item's data.
    A[index].d->m = -4294967296.0;
    A[index].d->WT_obj = -4294967296.0;
    A[index].d->MB_obj = -4294967296.0;
    A[index].d->i = 0;
    A[index].d->j = 0;

    // Swap this element with item at end of heap.
    tuple_o* swap = A[index].d;
    A[index].d = A[small_val].d;
    A[small_val].d = swap;
    A[index].d->k = index;      // Note the change in locations.
    A[small_val].d->k = small_val;

    // Note change in heap size.
    heaplimit--;

    // Downsift moved element to new location; O(log k).
    downsift(index);

    // Shrink by factor of 2 if necessary.
    if ((heaplimit / 2 > heapmin) && (heaplimit == (arraysize / 2))) shrink();
  }
}

bool maxheap::heapIsEmpty() { return isempty; }
int maxheap::heapSize() { return heaplimit - 1; }

void maxheap::printHeap()
{
  for (int i = 1; i < heaplimit; i++)
  {
#if defined(WT)
    std::cout << A[i].d << "\t[" << A[i].d->k << "]\tdQ = " << A[i].d->m
              << "\tWT_obj = " << A[i].d->WT_obj << "\t" << A[i].d->i << " -> "
              << A[i].d->j << "\n";
#elif defined(MB)
    std::cout << A[i].d << "\t[" << A[i].d->k << "]\tdQ = " << A[i].d->m
               << "\tMB_obj = " << A[i].d->MB_obj << "\t" << A[i].d->i << " -> "
               << A[i].d->j << "\n";
#else
    std::cout << A[i].d << "\t[" << A[i].d->k << "]\tdQ = " << A[i].d->m << "\t"
              << A[i].d->i << " -> " << A[i].d->j << "\n";
#endif
  }
}

void maxheap::printHeapTop10()
{
  int limit = heaplimit > 10 ? 11 : heaplimit;

  for (int i = 1; i < limit; i++)
  {
#if defined(WT)
    std::cout << A[i].d << "\t[" << A[i].d->k << "]\tdQ = " << A[i].d->m
              << "\tWT_obj = " << A[i].d->WT_obj << "\t" << A[i].d->i << " -> "
              << A[i].d->j << "\n";
#elif defined(MB)
    std::cout << A[i].d << "\t[" << A[i].d->k << "]\tdQ = " << A[i].d->m
              << "\tMB_obj = " << A[i].d->MB_obj << "\t" << A[i].d->i << " -> "
              << A[i].d->j << "\n";
#else
    std::cout << A[i].d << "\t[" << A[i].d->k << "]\tdQ = " << A[i].d->m << "\t"
              << A[i].d->i - 1 << " -> " << A[i].d->j - 1 << "\n";
#endif
  }

  // Changed slightly from previous version; suggested by Janne Aukia
  // (jaukia@cs.hut.fi) on 12 Oct 2006.
}

#endif

}
