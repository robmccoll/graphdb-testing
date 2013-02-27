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
/*! \file test_hachar.cpp

    \author Greg Mackey (gemacke@sandia.gov)

    \brief Tests the functionality of the hachar hash table.

    \date 7/14/2010
*/
/****************************************************************************/

#include <mtgl/hachar.hpp>
#include <mtgl/dynamic_array.hpp>
#include <mtgl/util.hpp>

using namespace mtgl;

typedef hachar<int, int> hash_table_type;

#define MAX_LOOP 10000000
#define HASH_SIZE 10000000

template<typename HT>
class table_visitor {
public:
  void operator()(typename HT::value_type* i) { i->second += i->first; }
};

template<typename V>
class insert_visitor {
public:
  insert_visitor(dynamic_array<V>& a) : arr(a) {}

  void operator()(V& val) { arr[val] = val; }

private:
  dynamic_array<V>& arr;
};

int main()
{
  int data;
  bool truth;
  table_visitor<hash_table_type> tv1;
  hash_mt_incr<int> uv1;

  dynamic_array<int> arr(MAX_LOOP);
  insert_visitor<int> iv1(arr);

  mt_timer timer;

  #pragma mta fence
  timer.start();
  hash_table_type xht(HASH_SIZE);
  #pragma mta fence
  timer.stop();

  printf("                 MAX_LOOP: %9d\n\n", MAX_LOOP);

  printf("           Initialization: %9.6lf\n", timer.getElapsedSeconds());

  #pragma mta fence
  timer.start();
  #pragma mta assert parallel
  #pragma mta block dynamic schedule
  for (int i = 0; i < MAX_LOOP; ++i)
  {
//    if (!xht.insert(i, i).second)  printf("Couldn't insert %d.\n", i);
//    xht.insert(i, i);
//    hash_table_type::value_type* ptr = xht.insert(i, i).first;
//    arr[ptr->second] = ptr->second;
    int v = xht.insert(i, i).first->second;
    arr[v] = v;
//    xht.insert(i, i, iv1);
  }
  #pragma mta fence
  timer.stop();

  printf("                Insertion: %9.6lf\n", timer.getElapsedSeconds());

/*
  printf("\n");
  xht.print_array_chunk_structure();
  printf("\n");
  xht.print_stats();
  printf("\n");
  printf("             Num Elements: %9lu\n", xht.size());
  printf("\n");
*/

  #pragma mta fence
  timer.start();
  hash_table_type xht2(xht);
  #pragma mta fence
  timer.stop();

  printf("         Copy constructor: %9.6lf\n", timer.getElapsedSeconds());

  #pragma mta fence
  timer.start();
  xht2 = xht;
  #pragma mta fence
  timer.stop();

  printf("      Assignment operator: %9.6lf\n", timer.getElapsedSeconds());

//  if (xht2.lookup(20, data))  printf("Found: 20 -> %d\n", data);
//  else  printf("Value 20 not found.\n");

  #pragma mta fence
  timer.start();
  xht2.visit(tv1);
  #pragma mta fence
  timer.stop();

//  if (xht2.lookup(20, data))  printf("Found: 20 -> %d\n", data);
//  else  printf("Value 20 not found.\n");

  printf("                    Visit: %9.6lf\n", timer.getElapsedSeconds());

  #pragma mta fence
  timer.start();
  xht2.clear();
  #pragma mta fence
  timer.stop();

  printf("                    Clear: %9.6lf\n", timer.getElapsedSeconds());

  #pragma mta fence
  timer.start();
  #pragma mta assert parallel
  #pragma mta block dynamic schedule
  for (int i = 0; i < MAX_LOOP; ++i)
  {
    xht2.insert(i, i, iv1);
  }
  #pragma mta fence
  timer.stop();

  printf("        Visitor Insertion: %9.6lf\n", timer.getElapsedSeconds());

/*
  printf("\narr:");
  for (int i = 0; i < MAX_LOOP; ++i) printf("%4d", arr[i]);
  printf("\n\n");
*/

/*
  #pragma mta fence
  timer.start();
  xht2.resize(10000);
  #pragma mta fence
  timer.stop();

  printf("             Empty Resize: %9.6lf\n", timer.getElapsedSeconds());
*/

  #pragma mta fence
  timer.start();
  #pragma mta block dynamic schedule
  #pragma mta assert nodep
  for (int i = 0; i < MAX_LOOP; ++i)
  {
//    if (!xht.member(i))  printf("Didn't find %d.\n", i);
    truth = xht.member(i);
  }
  #pragma mta fence
  timer.stop();

  printf("          Membership test: %9.6lf        %d\n",
         timer.getElapsedSeconds(), truth);

  #pragma mta fence
  timer.start();
  #pragma mta block dynamic schedule
  #pragma mta assert nodep
  for (int i = 0; i < MAX_LOOP; ++i)
  {
//    if (!xht.lookup(i, data))  printf("Didn't find %d.\n", i);
    data = xht.lookup(i, data);
  }
  #pragma mta fence
  timer.stop();

  printf("                   Lookup: %9.6lf        %d\n",
         timer.getElapsedSeconds(), data);

/*
  #pragma mta fence
  timer.start();
  #pragma mta assert parallel
  #pragma mta block dynamic schedule
  for (int i = 0; i < MAX_LOOP; ++i)
  {
    xht[i] = i + 1;
  }
  #pragma mta fence
  timer.stop();

  printf("             Write via []: %9.6lf        %d\n",
         timer.getElapsedSeconds(), data);
*/

  #pragma mta fence
  timer.start();
  int sum = 0;
  #pragma mta assert parallel
  #pragma mta assert nodep
  #pragma mta block dynamic schedule
  for (int i = 0; i < MAX_LOOP; ++i)
  {
    sum += xht[i];
  }
  #pragma mta fence
  timer.stop();

  printf("              Read via []: %9.6lf        %d\n",
         timer.getElapsedSeconds(), sum);

/*
  #pragma mta fence
  timer.start();
  #pragma mta assert parallel
  #pragma mta block dynamic schedule
  for (int i = 0; i < MAX_LOOP; ++i)
  {
    xht[i] = xht[i] + 4;
  }
  #pragma mta fence
  timer.stop();

  printf("         Usage and lookup: %9.6lf\n", timer.getElapsedSeconds());
*/

  #pragma mta fence
  timer.start();
  #pragma mta block dynamic schedule
  #pragma mta assert nodep
  for (int i = 0; i < MAX_LOOP; ++i)
  {
//    if (!xht.update(i, i * 2))  printf("Couldn't update %d.\n", i);
    xht.update(i, i * 2);
  }
  #pragma mta fence
  timer.stop();

  printf("                   Update: %9.6lf\n", timer.getElapsedSeconds());

  #pragma mta fence
  timer.start();
  #pragma mta block dynamic schedule
  #pragma mta assert nodep
  for (int i = 0; i < MAX_LOOP; ++i)
  {
//    if (!xht.update(i, i * 2))  printf("Couldn't update %d.\n", i);
    xht.update(i, 2, uv1);
  }
  #pragma mta fence
  timer.stop();

  printf("           Update-mt_incr: %9.6lf\n", timer.getElapsedSeconds());

  #pragma mta fence
  timer.start();
  #pragma mta block dynamic schedule
  #pragma mta assert parallel
  for (int i = 0; i < MAX_LOOP; ++i)
  {
    xht.update_insert(i, i * 3);
  }
  #pragma mta fence
  timer.stop();

  printf("            Update-Insert: %9.6lf\n", timer.getElapsedSeconds());

  #pragma mta fence
  timer.start();
  #pragma mta block dynamic schedule
  #pragma mta assert parallel
  for (int i = 0; i < MAX_LOOP; ++i)
  {
    xht.update_insert(i, 3, uv1);
  }
  #pragma mta fence
  timer.stop();

  printf("    Update-Insert-mt_incr: %9.6lf\n", timer.getElapsedSeconds());

  #pragma mta fence
  timer.start();
  assign_contiguous_ids(xht);
  #pragma mta fence
  timer.stop();

  printf("    Assign contiguous ids: %9.6lf\n", timer.getElapsedSeconds());

/*
  #pragma mta fence
  timer.start();
  xht.resize(MAX_LOOP);
  #pragma mta fence
  timer.stop();

  printf("              Full Resize: %9.6lf\n", timer.getElapsedSeconds());

//  if (xht.lookup(20, data))  printf("Found: 20 -> %d\n", data);

  #pragma mta fence
  timer.start();
  #pragma mta block dynamic schedule
  #pragma mta assert parallel
  for (int i = 0; i < MAX_LOOP; i += 2)
  {
//    if (!xht.erase(i))  printf("Didn't find %d.\n", i);
    xht.erase(i);
  }
  #pragma mta fence
  timer.stop();

  printf("                    Erase: %9.6lf\n", timer.getElapsedSeconds());
*/

//  if (xht.lookup(20, data))  printf("Found: 20 -> %d\n", data);
//  else  printf("Value 20 not found.\n");
//  if (xht.lookup(139, data))  printf("Found: 139 -> %d\n", data);
//  else  printf("Value 139 not found.\n");

  printf("\n");
  printf("             Num Elements: %9lu\n", xht.size());

  return 0;
}
