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
/*! \file test_xmt_hash_table.cpp

    \author Greg Mackey (gemacke@sandia.gov)

    \brief Tests the functionality of the xmt_hash_table.

    \date 3/5/2010
*/
/****************************************************************************/

#include <mtgl/xmt_hash_table.hpp>
#include <mtgl/dynamic_array.hpp>
#include <mtgl/util.hpp>
#include <mtgl/random.hpp>

using namespace mtgl;

typedef xmt_hash_table<int, int> hash_table_type;

#define MAX_LOOP 10000000

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

int main(int argc, const char* argv[])
{
  unsigned long num = MAX_LOOP;
  int cur_arg_index = 1;
  bool print = false;

  while (cur_arg_index < argc)
  {
    if (strcmp(argv[cur_arg_index], "-num") == 0)
    {
      ++cur_arg_index;
      num = strtoul(argv[cur_arg_index], NULL, 10);
    }
    else if (strcmp(argv[cur_arg_index], "-print") == 0)
    {
      print = true;
    }
    else
    {
      printf("Unrecognized option: %s\n", argv[cur_arg_index]);
    }

    ++cur_arg_index;
  }

  int data = 0;
  bool truth = false;
  table_visitor<hash_table_type> tv1;
  hash_mt_incr<int> uv1;

  dynamic_array<int> arr(num);
  insert_visitor<int> iv1(arr);

  mt_timer timer;

  #pragma mta fence
  timer.start();
  hash_table_type xht(2 * num);
  #pragma mta fence
  timer.stop();

  printf("                      num: %9lu\n", num);
  printf("               capacity(): %9lu\n\n", xht.capacity());
  printf("           Initialization: %9.6lf\n", timer.getElapsedSeconds());

  #pragma mta fence
  timer.start();
  #pragma mta assert parallel
  #pragma mta block schedule
  for (unsigned long i = 0; i < num; ++i)
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
  for (unsigned long i = 0; i < num; ++i)
  {
    xht2.insert(i, i, iv1);
  }
  #pragma mta fence
  timer.stop();

  printf("        Visitor Insertion: %9.6lf\n", timer.getElapsedSeconds());

/*
  printf("\narr:");
  for (unsigned long i = 0; i < num; ++i) printf("%4d", arr[i]);
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
  #pragma mta block schedule
  for (unsigned long i = 0; i < num; ++i)
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
  #pragma mta block schedule
  for (unsigned long i = 0; i < num; ++i)
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
  #pragma mta block schedule
  for (unsigned long i = 0; i < num; ++i)
  {
    xht[i] = i + 1;
  }
  #pragma mta fence
  timer.stop();

  printf("             Write via []: %9.6lf        %d\n",
         timer.getElapsedSeconds(), data);

  #pragma mta fence
  timer.start();
  int sum = 0;
  #pragma mta assert parallel
  #pragma mta block schedule
  for (unsigned long i = 0; i < num; ++i)
  {
    sum += xht[i];
  }
  #pragma mta fence
  timer.stop();

  printf("              Read via []: %9.6lf        %d\n",
         timer.getElapsedSeconds(), sum);

  #pragma mta fence
  timer.start();
  #pragma mta assert parallel
  #pragma mta block schedule
  for (unsigned long i = 0; i < num; ++i)
  {
    xht[i] = xht[i] + 4;
  }
  #pragma mta fence
  timer.stop();

  printf("         Usage and lookup: %9.6lf\n", timer.getElapsedSeconds());
*/

  #pragma mta fence
  timer.start();
  #pragma mta block schedule
  #pragma mta assert nodep
  for (unsigned long i = 0; i < num; ++i)
  {
//    if (!xht.update(i, i * 2))  printf("Couldn't update %d.\n", i);
    xht.update(i, i * 2);
  }
  #pragma mta fence
  timer.stop();

  printf("                   Update: %9.6lf\n", timer.getElapsedSeconds());

  #pragma mta fence
  timer.start();
  #pragma mta block schedule
  #pragma mta assert nodep
  for (unsigned long i = 0; i < num; ++i)
  {
//    if (!xht.update(i, i * 2))  printf("Couldn't update %d.\n", i);
    xht.update(i, 2, uv1);
  }
  #pragma mta fence
  timer.stop();

  printf("           Update-mt_incr: %9.6lf\n", timer.getElapsedSeconds());

  #pragma mta fence
  timer.start();
  #pragma mta block schedule
  #pragma mta assert parallel
  for (unsigned long i = 0; i < num; ++i)
  {
    xht.update_insert(i, i * 3);
  }
  #pragma mta fence
  timer.stop();

  printf("            Update-Insert: %9.6lf\n", timer.getElapsedSeconds());

  #pragma mta fence
  timer.start();
  #pragma mta block schedule
  #pragma mta assert parallel
  for (unsigned long i = 0; i < num; ++i)
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

  #pragma mta fence
  timer.start();
  xht.resize(4 * num);
  #pragma mta fence
  timer.stop();

  printf("              Full Resize: %9.6lf\n", timer.getElapsedSeconds());

//  if (xht.lookup(20, data))  printf("Found: 20 -> %d\n", data);

  #pragma mta fence
  timer.start();
  #pragma mta block schedule
  #pragma mta assert parallel
  for (unsigned long i = 0; i < num; i++)
  {
//    if (!xht.erase(i))  printf("Didn't find %d.\n", i);
    xht.erase(i);
  }
  #pragma mta fence
  timer.stop();

  printf("                    Erase: %9.6lf\n", timer.getElapsedSeconds());

  unsigned long table_size = xht.size();
  printf("\nChecking to see that we erased everything.\n");
  if (table_size != 0)
  {
    printf("\tERROR with erasing, expected 0 elements but found %lu in "
           "the table.\n", table_size);
  }
  else
  {
    printf("\tSUCCESS!\n");
  }

//  if (xht.lookup(20, data))  printf("Found: 20 -> %d\n", data);
//  else  printf("Value 20 not found.\n");
//  if (xht.lookup(139, data))  printf("Found: 139 -> %d\n", data);
//  else  printf("Value 139 not found.\n");

  printf("\n");
  printf("             Num Elements: %9lu\n", table_size);
  printf("               capacity(): %9lu\n\n", xht.capacity());

  if (print) xht.print();

  timer.start();
  #pragma mta block schedule
  #pragma mta assert parallel
  for (unsigned long i = 0; i < num; i++)
  {
//    if (!xht.erase(i))  printf("Didn't find %d.\n", i);
    xht.insert_reuse(i, i);
  }
  #pragma mta fence
  timer.stop();

  if (print) xht.print();
  printf("              Insert reuse: %9.6lf\n", timer.getElapsedSeconds());

  // Checking summation of values
  unsigned long total = 0;
  unsigned long num_streams = 1;
  unsigned long stream_id = 0;
  #pragma mta for all streams stream_id of num_streams
  {
    hash_table_type::thread_iterator begin_iter = xht.begin(stream_id,
                                                            num_streams);
    hash_table_type::thread_iterator end_iter = xht.end(stream_id, num_streams);

    unsigned long local_total = 0;

    for ( ; begin_iter != end_iter; ++begin_iter)
    {
      local_total += begin_iter->second;
    }

    mt_incr(total, local_total);
  }
  #pragma mta fence

  unsigned long expected_total = num * (num - 1) / 2;

  printf("\n");
  printf("Checking to see that the expected sum of the values equals"
         " the actual sum\n");
  if (expected_total != total)
  {
    printf("\tERROR: actual sum %lu does not match expected sum %lu\n\n",
            total, expected_total);
  }
  else
  {
    printf("\tSUCCESS!\n\n");
  }

  // This next test we erase everything again and then add random integers
  // instead of consecutive integers.

  #pragma mta fence
  timer.start();
  #pragma mta block schedule
  #pragma mta assert parallel
  for (unsigned long i = 0; i < num; i++)
  {
    xht.erase(i);
  }
  #pragma mta fence
  timer.stop();

  printf("              Erase Again: %9.6lf\n\n", timer.getElapsedSeconds());
  table_size = xht.size();
  printf("Checking to see that we erased everything.\n");
  if (table_size != 0)
  {
    printf("\tERROR with erasing, expected %d elements but found %lu in "
           "the table.\n", 0, table_size);
  }
  else
  {
    printf("\tSUCCESS!\n");
  }

  printf("\n");
  printf("             Num Elements: %9lu\n", xht.size());
  printf("               capacity(): %9lu\n\n", xht.capacity());

  if (print) xht.print();

  return 0;
}
