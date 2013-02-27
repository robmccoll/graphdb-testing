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
/*! \file vektor.h

    \brief This class is only here because it is needed by wcnm().  It
           should not be used for any other purpose.
*/
/****************************************************************************/

#if !defined(vektor_INCLUDED)
#define vektor_INCLUDED

#include <iostream>

#if !defined(vektor_INCLUDED)
#define vektor_INCLUDED
#include <mtgl/maxheap.h>
#endif

namespace mtgl {

#if !defined(DPAIR_INCLUDED)
#define DPAIR_INCLUDED
class dpair {
public:
  int x;
  double y;
  dpair* next;

  dpair() : x(0), y(0.0), next(NULL) {}
  ~dpair() {}
};
#endif

struct dppair {
  dpair* head;
  dpair* tail;
};

class element {
public:
  int key;             ///< Binary-tree key.
  double stored;       ///< Additional stored value (associated with key).
  tuple_o* heap_ptr;   ///< Pointer to element's location in vektor max-heap.

  bool color;          ///< F: BLACK, T: RED.

  element* parent;     ///< Pointer to parent node.
  element* left;       ///< Pointer for left subtree.
  element* right;      ///< Pointer for right subtree.

  element() : key(0), stored(-4294967296.0), color(false),
              parent(NULL), left(NULL), right(NULL) {}
  ~element() {}
};

/*!  This vector implementation is a pair of linked data structures: a
     red-black balanced binary tree data structure and a maximum heap. This
     pair allows us to find a stored element in time O(log n), find the
     maximum element in time O(1), update the maximum element in time
     O(log n), delete an element in time O(log n), and insert an element in
     time O(log n). These complexities allow a much faster implementation of
     the fastCommunity algorithm. If we dispense with the max-heap, then some
     operations related to updating the maximum stored value can take up to
     O(n), which is potentially very slow.

     Both the red-black balanced binary tree and the max-heap implementations
     are custom-jobs. Note that the key=0 is assumed to be a special value,
     and thus you cannot insert such an item. Beware of this limitation.
*/
class vektor {
private:
  int my_id;           ///< Will be a node id in CNM.
  element* root;       ///< Binary tree root.
  element* leaf;       ///< All leaf nodes.
  maxheap* heap;       ///< Max-heap of elements in vektor.
  int support;         ///< Number of nodes in the tree.

  /// Left-rotation operator.
  void rotateLeft(element* x);

  /// Right-rotation operator.
  void rotateRight(element* y);

  /// House-keeping after insertion.
  void insertCleanup(element* z);

  /// House-keeping after deletion.
  void deleteCleanup(element* x);

  /// Internal recursive cons'ing function.
  dppair* consSubtree(element* z);

  dpair* returnSubtreeAsList(element* z, dpair* head);

  /// Display the subtree rooted at z.
  void printSubTree(element* z);

  /// Delete subtree rooted at z.
  void deleteSubTree(element* z);

  /// Returns minimum of subtree rooted at z.
  element* returnMinKey(element* z);

  /// Returns successor of z's key.
  element* returnSuccessor(element* z);

public:
  /// Default constructor.
  vektor(int id, int size);

  /// Default destructor.
  ~vektor();

  /// Returns T if searchKey found, and points foundNode at the
  /// corresponding node.
  element* findItem(const int searchKey);

  /// Insert a new key with stored value.
  void insertItem(int newKey, double newStored);

  /// Selete a node with given key.
  void deleteItem(int killKey);

  /// Delete the entire tree.
  void deleteTree();

  /// Return the tree as a list of dpairs.
  dpair* returnTreeAsList();

  /// Return the tree as a list of dpairs.
  dpair* returnTreeAsList2();

  /// Returns the maximum key in the tree.
  tuple_o returnMaxKey();

  /// Returns a tuple of the maximum (key, .stored).
  tuple_o returnMaxStored();

  /// Returns number of items in tree.
  int returnNodecount();

  /// Displays tree (in-order traversal).
  void printTree();

  /// Displays heap.
  void printHeap();

  int returnArraysize();
  int returnHeaplimit();
};


vektor::vektor(int id, int size) : my_id(id)
{
  root = new element;
  leaf = new element;
  heap = new maxheap(size);

  leaf->parent = root;

  root->left = leaf;
  root->right = leaf;
  support = 0;
}

vektor::~vektor()
{
  if (root != NULL) deleteSubTree(root);

  support = 0;
  delete leaf;
  root = NULL;
  leaf = NULL;
  delete heap;
  heap = NULL;
}

void vektor::deleteTree()
{
  if (root != NULL) deleteSubTree(root);
}

void vektor::deleteSubTree(element* z)
{
  if (z->left != leaf) deleteSubTree(z->left);
  if (z->right != leaf) deleteSubTree(z->right);

  delete z;
  z = NULL;
}


/// Public search function - if there exists a element in the three with
/// key=searchKey, it returns TRUE and foundNode is set to point to the found
/// node; otherwise, it sets foundNode=NULL and returns FALSE.
element* vektor::findItem(const int searchKey)
{
  element* current = root;

  // Empty tree; bail out.
  if (current->key == 0) return NULL;

  while (current != leaf)
  {
    // Left or right?
    if (searchKey < current->key)
    {
      // Try moving down-left.
      if (current->left != leaf)
      {
        current = current->left;
      }
      else
      {
        // Failure; bail out.
        return NULL;
      }
    }
    else
    {
      // Left or right?
      if (searchKey > current->key)
      {
        // Try moving down-left.
        if (current->right != leaf)
        {
          current = current->right;
        }
        else
        {
          // Failure; bail out.
          return NULL;
        }
      }
      else
      {
        // Found (searchKey==current->key).
        return current;
      }
    }
  }

  return NULL;
}

dpair* vektor::returnTreeAsList()   // pre-order traversal
{
  dpair* head = new dpair;
  head->x = root->key;
  head->y = root->stored;

  dpair* tail = head;

  if (root->left != leaf) tail = returnSubtreeAsList(root->left, tail);
  if (root->right != leaf) tail = returnSubtreeAsList(root->right, tail);

  return head->x == 0 ? NULL : head;
}

dpair* vektor::returnSubtreeAsList(element* z, dpair* head)
{
  dpair* newnode = new dpair;
  newnode->x = z->key;
  newnode->y = z->stored;

  head->next = newnode;

  dpair* tail = newnode;

  if (z->left != leaf) tail = returnSubtreeAsList(z->left, tail);
  if (z->right != leaf) tail = returnSubtreeAsList(z->right, tail);

  return tail;
}

tuple_o vektor::returnMaxStored()
{
  return heap->returnMaximum();
}

tuple_o vektor::returnMaxKey()
{
  element* current = root;

  // Search to bottom-right corner of tree.
  while (current->right != leaf) current = current->right;

  // Store the data found.
  tuple_o themax;
  themax.m = current->stored;
  themax.i = current->key;
  themax.j = current->key;

  // Return that data.
  return themax;
}

element* vektor::returnMinKey(element* z)
{
  element* current = z;

  // Search to bottom-right corner of tree.
  while (current->left != leaf) current = current->left;

  // Return pointer to the minimum.
  return current;
}

element* vektor::returnSuccessor(element* z)
{
  element* w = z;

  // If right-subtree exists, return min of it.
  if (w->right != leaf) return returnMinKey(w->right);

  // Else search up in tree.
  element* current = w->parent;
  while ((current != NULL) && (w == current->right))
  {
    // Move up in tree until find a non-right-child.
    w = current;
    current = current->parent;
  }

  return current;
}

int vektor::returnNodecount() { return support; }
int vektor::returnArraysize() { return heap->returnArraysize(); }
int vektor::returnHeaplimit() { return heap->returnHeaplimit(); }

void vektor::insertItem(int newKey, double newStored)
{
  // First we check to see if newKey is already present in the tree. If so,
  // we simply set .stored += newStored. If not, we must find where to insert
  // the key.

  // Find newKey in tree; return pointer to it O(log k).
  element* current = findItem(newKey);
  if (current != NULL)
  {
    // JWB: This used to be incr.
    current->stored = newStored;

    // Update corresponding element in heap + reheapify; O(log k).
    heap->updateItem(current->heap_ptr, current->stored);
  }
  else
  {
    // Didn't find it, so need to create it.

    tuple_o newitem;
    newitem.m = newStored;
    newitem.i = my_id;
    newitem.j = newKey;

    element* newNode = new element;  // Element for the vektor.
    newNode->key = newKey;           // Store newKey.
    newNode->stored = newStored;     // Store newStored.
    newNode->color = true;           // New nodes are always RED.
    newNode->parent = NULL;          // New node initially has no parent.
    newNode->left = leaf;            // Left leaf.
    newNode->right = leaf;           // Right leaf.

    // Add new item to the vektor heap.
    newNode->heap_ptr = heap->insertItem(newitem);

    support++;                       // Increment node count in vektor.

    // Must now search for where to insert newNode, i.e., find the correct
    // parent and set the parent and child to point to each other properly.
    current = root;

    if (current->key == 0)
    {
      // Insert as root.
      delete root;                   // Delete old root.
      root = newNode;                // Set root to newNode.
      leaf->parent = newNode;        // Set leaf's parent.
      current = leaf;                // Skip next loop.
    }

    // Search for insertion point.
    while (current != leaf)
    {
      // Left or right?
      if (newKey < current->key)
      {
        // Try moving down-left.
        if (current->left != leaf)
        {
          current = current->left;
        }
        else
        {
          // Found new parent.
          newNode->parent = current; // Set parent.
          current->left = newNode;   // Set child.
          current = leaf;            // Exit search.
        }
      }
      else
      {
        // Try moving down-right.
        if (current->right != leaf)
        {
          current = current->right;
        }
        else
        {
          // Found new parent.
          newNode->parent = current; // Set parent.
          current->right = newNode;  // Set child.
          current = leaf;            // Exit search.
        }
      }
    }

    // Now do the house-keeping necessary to preserve the red-black
    // properties and maintain balance.
    insertCleanup(newNode);
  }
}

void vektor::insertCleanup(element* z)
{
  // Fix now if z is root.
  if (z->parent == NULL)
  {
    z->color = false;
    return;
  }

  // While z is not root and z's parent is RED.
  while (z->parent != NULL && z->parent->color)
  {
    if (z->parent == z->parent->parent->left)
    {
      // z's parent is LEFT-CHILD.
      element* temp = z->parent->parent->right;   // Grab z's uncle.

      if (temp->color)
      {
        z->parent->color = false;        // Color z's parent BLACK.    (Case 1)
        temp->color = false;             // Color z's uncle BLACK.     (Case 1)
        z->parent->parent->color = true; // Color z's grandparent RED. (Case 1)
        z = z->parent->parent;           // Set z = z's grandparent.   (Case 1)
      }
      else
      {
        // z is RIGHT-CHILD.
        if (z == z->parent->right)
        {
          z = z->parent;                 // Set z = z's parent.        (Case 2)
          rotateLeft(z);                 // Perform left-rotation.     (Case 2)
        }

        z->parent->color = false;        // Color z's parent BLACK.    (Case 3)
        z->parent->parent->color = true; // Color z's grandparent RED. (Case 3)
        rotateRight(z->parent->parent);  // Perform right-rotation.    (Case 3)
      }
    }
    else
    {
      // z's parent is RIGHT-CHILD.
      element* temp = z->parent->parent->left;    // Grab z's uncle.

      if (temp->color)
      {
        z->parent->color = false;        // Color z's parent BLACK.    (Case 1)
        temp->color = false;             // Color z's uncle BLACK.     (Case 1)
        z->parent->parent->color = true; // Color z's grandparent RED. (Case 1)
        z = z->parent->parent;           // Set z = z's grandparent.   (Case 1)
      }
      else
      {
        // z is LEFT-CHILD.
        if (z == z->parent->left)
        {
          z = z->parent;                 // Set z = z's parent.        (Case 2)
          rotateRight(z);                // Perform right-rotation.    (Case 2)
        }

        z->parent->color = false;        // Color z's parent BLACK.    (Case 3)
        z->parent->parent->color = true; // Color z's grandparent RED. (Case 3)
        rotateLeft(z->parent->parent);   // Perform left-rotation.     (Case 3)
      }
    }
  }

  // Color the root BLACK.
  root->color = false;
}

void vektor::deleteItem(int killKey)
{
  element* z = findItem(killKey);

  // Item not present; bail out.
  if (z == NULL) return;

  if (z != NULL)
  {
    // Get old maximum in O(1).
    tuple_o newmax = heap->returnMaximum();

    // Delete item in the max-heap O(log k).
    heap->deleteItem(z->heap_ptr);
  }

  // Attempt to delete the root.
  if (support == 1)
  {
    // Restore root node to default state.
    root->key = 0;
    root->stored = -4294967296.0;
    root->color = false;
    root->parent = NULL;
    root->left = leaf;
    root->right = leaf;
    root->heap_ptr = NULL;

    // Set support to zero.
    support--;

    // Exit - no more work to do.
    return;
  }

  if (z != NULL)
  {
    element* x;
    element* y;

    // Decrement node count.
    support--;

    // Case of less than two children
    if ((z->left == leaf) || (z->right == leaf))
    {
      // Set y to be z.
      y = z;
    }
    else
    {
      // Set y to be z's key-successor.
      y = returnSuccessor(z);
    }

    if (y->left != leaf)
    {
      // Pick y's one child (left-child).
      x = y->left;
    }
    else
    {
      // Right-child.
      x = y->right;
    }

    // Make y's child's parent be y's parent.
    x->parent = y->parent;

    if (y->parent == NULL)
    {
      // If y is the root, x is now root.
      root = x;
    }
    else
    {
      // Decide y's relationship with y's parent.
      if (y == y->parent->left)
      {
        // Replace x as y's parent's left child.
        y->parent->left = x;
      }
      else
      {
        // Replace x as y's parent's left child.
        y->parent->right = x;
      }
    }

    // Insert y into z's spot.
    if (y != z)
    {
      // Copy y data into z.
      z->key = y->key;
      z->stored = y->stored;
      z->heap_ptr = y->heap_ptr;
    }

    // Do house-keeping to maintain balance.
    if (y->color == false) deleteCleanup(x);

    // Deallocate y.
    delete y;

    // Point y to NULL for safety.
    y = NULL;
  }
}

void vektor::deleteCleanup(element* x)
{
  // Until x is the root, or x is RED.
  while ((x != root) && (x->color == false))
  {
    if (x == x->parent->left)
    {
      // Branch on x being a LEFT-CHILD.
      element* w = x->parent->right; // Grab x's sibling.

      if (w->color == true)          // If x's sibling is RED.
      {
        w->color = false;            // Color w BLACK.                 (Case 1)
        x->parent->color = true;     // Color x's parent RED.          (Case 1)
        rotateLeft(x->parent);       // Left rotation on x's parent.   (Case 1)
        w = x->parent->right;        // Make w be x's right sibling.   (Case 1)
      }

      if ((w->left->color == false) && (w->right->color == false))
      {
        w->color = true;             // Color w RED.                   (Case 2)
        x = x->parent;               // Examine x's parent.            (Case 2)
      }
      else
      {
        if (w->right->color == false)
        {
          w->left->color = false;    // Color w's left child BLACK.    (Case 3)
          w->color = true;           // Color w RED.                   (Case 3)
          element* t = x->parent;    // Store x's parent.
          rotateRight(w);            // Right rotation on w.           (Case 3)
          x->parent = t;             // Restore x's parent.
          w = x->parent->right;      // Make w be x's right sibling.   (Case 3)
        }

        w->color = x->parent->color; // Make w's color = x's parent's. (Case 4)
        x->parent->color = false;    // Color x's parent BLACK.        (Case 4)
        w->right->color = false;     // Color w's right child BLACK.   (Case 4)
        rotateLeft(x->parent);       // Left rotation on x's parent.   (Case 4)
        x = root;                    // Finished work. bail out.       (Case 4)
      }
    }
    else
    {
      // x is RIGHT-CHILD.
      element* w = x->parent->left;  // Grab x's sibling.

      if (w->color == true)          // If x's sibling is RED.
      {
        w->color = false;            // Color w BLACK.                 (Case 1)
        x->parent->color = true;     // Color x's parent RED.          (Case 1)
        rotateRight(x->parent);      // Right rotation on x's parent.  (Case 1)
        w = x->parent->left;         // Make w be x's left sibling.    (Case 1)
      }

      if ((w->right->color == false) && (w->left->color == false))
      {
        w->color = true;             // Color w RED.                   (Case 2)
        x = x->parent;               // Examine x's parent.            (Case 2)
      }
      else
      {
        if (w->left->color == false)
        {
          w->right->color = false;   // Color w's right child BLACK.   (Case 3)
          w->color = true;           // Color w RED.                   (Case 3)
          element* t = x->parent;    // Store x's parent.
          rotateLeft(w);             // Left rotation on w.            (Case 3)
          x->parent = t;             // Restore x's parent.
          w = x->parent->left;       // Make w be x's left sibling.    (Case 3)
        }

        w->color = x->parent->color; // Make w's color = x's parent's. (Case 4)
        x->parent->color = false;    // Color x's parent BLACK.        (Case 4)
        w->left->color = false;      // Color w's left child BLACK.    (Case 4)
        rotateRight(x->parent);      // Right rotation on x's parent.  (Case 4)
        x = root;                    // x is now the root.             (Case 4)
      }
    }
  }

  // Color x (the root) BLACK.
  x->color = false;
}

void vektor::rotateLeft(element* x)
{
  // Do pointer-swapping operations for left-rotation.
  element* y = x->right;         // Grab right child.
  x->right = y->left;            // Make x's RIGHT-CHILD be y's LEFT-CHILD.
  y->left->parent = x;           // Make x be y's LEFT-CHILD's parent.
  y->parent = x->parent;         // Make y's new parent be x's old parent.

  if (x->parent == NULL)
  {
    // If x was root, make y root.
    root = y;
  }
  else
  {
    // If x is LEFT-CHILD, make y be x's parent's.
    if (x == x->parent->left)
    {
      // Left-child.
      x->parent->left = y;
    }
    else
    {
      // Right-child.
      x->parent->right = y;
    }
  }

  // Make x be y's LEFT-CHILD.
  y->left = x;

  // Make y be x's parent.
  x->parent = y;
}

void vektor::rotateRight(element* y)
{
  // Do pointer-swapping operations for right-rotation.
  element* x = y->left;          // Grab left child
  y->left = x->right;            // Replace left child yith x's right subtree.
  x->right->parent = y;          // Replace y as x's right subtree's parent.
  x->parent = y->parent;         // Make x's new parent be y's old parent.

  if (y->parent == NULL)
  {
    // If y was root, make x root.
    root = x;
  }
  else
  {
    // If y is RIGHT-CHILD, make x be y's parent's.
    if (y == y->parent->right)
    {
      // Right-child.
      y->parent->right = x;
    }
    else
    {
      // Left-child.
      y->parent->left = x;
    }
  }

  // Make y be x's RIGHT-CHILD.
  x->right = y;

  // Make x be y's parent.
  y->parent = x;
}

void vektor::printTree()
{
  tuple_o max = heap->returnMaximum();
  std::cout << "\nTREE SIZE = " << support << std::endl;
  std::cout << "HEAP SIZE = " << heap->heapSize() << std::endl;
  std::cout << "MAXIMUM (" << max.j << " " << max.m << ")" << std::endl;
  std::cout << "# ";
  printSubTree(root);
}

void vektor::printHeap() { heap->printHeap(); }

void vektor::printSubTree(element* z)
{
  if (z == leaf) { return; }else
  {
    std::cout << "(" << z->key << " " << z->stored << " " << z->color << ") ["
         << z->heap_ptr << "]" << std::endl;

    std::cout << "L ";
    printSubTree(z->left);
    std::cout << std::endl;

    std::cout << "R ";
    printSubTree(z->right);
    std::cout << std::endl;
  }
}

}

#endif
