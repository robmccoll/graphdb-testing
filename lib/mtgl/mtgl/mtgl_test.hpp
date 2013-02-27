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
/*! \file mtgl_test.hpp

    \brief Provides a standard testing harness for use in test files.

    \author Greg Mackey (gemacke@sandia.gov)

    \date 4/13/2011

    Provides two functions that assist with testing.
    
    The function init_test() checks that the command-line arguments are
    correct and prints an error message if they weren't.  It also initializes
    the qthreads library if qthreads is being used.  There are two forms of
    the function that depend on whether or not the test takes user-defined
    arguments.  This function needs to be called before any parallel sections.
    
    The function create_test_graph() creates or reads in the graph specified
    by the command-line arguments.  There are two forms of the function that
    depend on whether or not to read graph weights.

    Note that currently only compressed_sparse_row_graph and
    stinger_graph_adapter support the mmap format.

    Here is an example of using the test harness with no user-defined
    arguments and no weights:

      int main(int argc, char *argv[])
      {
        init_test(argc, argv);

        compressed_sparse_row_graph<directedS> g;
        create_test_graph(g, argc, argv);
        ...
      }

    Here is an example of using the test harness with two user-defined
    arguments and weights:

      int main(int argc, char *argv[])
      {
        const int num_ts_args = 2;
        const char* ts_arg_names[num_ts_args] = { "arg1", "arg2" };
        const char* ts_arg_descs[num_ts_args] =
        {
          "The description for argument 1.",
          "The description for argument 2."
        };
        char** ts_argv;
        int ts_argc;
   
        init_test(argc, argv, num_ts_args, ts_arg_names,
                  ts_arg_descs, ts_argc, ts_argv);

        compressed_sparse_row_graph<directedS> g;
        dynamic_array<double> weights;
        create_test_graph(g, weights, argc, argv);
   
        int arg1 = atoi(ts_argv[0]);
        int arg2 = atoi(ts_argv[1]);
        ...
      }

*/
/****************************************************************************/

#ifndef MTGL_MTGL_TEST_HPP
#define MTGL_MTGL_TEST_HPP

#include <mtgl/mtgl_io.hpp>
#include <mtgl/generate_rmat_graph.hpp>
#include <mtgl/generate_erdos_renyi_graph.hpp>
#include <mtgl/dynamic_array.hpp>

#include <cstring>
#include <cctype>
#include <cstdlib>

#define MAX_DESC_LEN 2000
#define MAX_FORM_DESC_LEN 3 * MAX_DESC_LEN / 2

#ifdef USING_QTHREADS
#include <qthread.h>
#define NUM_FILE_ARGS 3
#define NUM_RMAT_ARGS 4
#define NUM_ER_ARGS 5
#elif defined(_OPENMP)
#include <omp.h>
#define NUM_FILE_ARGS 3
#define NUM_RMAT_ARGS 4
#define NUM_ER_ARGS 5
#else
#define NUM_FILE_ARGS 2
#define NUM_RMAT_ARGS 3
#define NUM_ER_ARGS 4
#endif

namespace mtgl {

namespace detail {

void print_test_description()
{
  fprintf(stderr,
    "DESCRIPTION\n"
    "    Specifying a generator uses the specified generator to generate a "
    "graph,\n"
#if defined(USING_QTHREADS) || defined (_OPENMP)
    "    and specifying a filename reads the graph from a file.  There will "
    "be\n"
    "    num_threads threads started.\n\n"
#else
    "    and specifying filename reads the graph from a file.\n\n"
#endif
    "GENERATORS\n"
    "    rmat p\n"
    "        Generates an R-MAT graph with 2^p vertices.\n\n"
    "    er n p\n"
    "        Generates an Erdos-Renyi graph with n vertices and edge "
    "probability\n"
    "        p.\n\n"
    "FILE FORMATS\n"
    "    '.dimacs'             Indicates the DIMACS format.\n"
    "    '.mtx'                Indicates the MatrixMarket format.\n"
    "    '.srcs' and '.dests'  Indicates the sources/dests binary format.\n"
    "    '.mmap'               Indicates a graph that has been mmapped "
    "into memory.\n");

  fprintf(stderr, "\n");
}

template <typename Graph, unsigned long mmap_status>
struct create_mmap_graph
{
  void operator()(Graph& g, char* argv[])
  {
    // Intialize with mmapped graph.
    read_mmap(g, argv[1]);
  }
};

template <typename Graph>
struct create_mmap_graph<Graph, MMAP_TYPE_NOT_DEFINED>
{
  void operator()(Graph& g, char* argv[])
  {
    fprintf(stderr, "This graph type doesn't support mmapping.\n");
    exit(1);
  }
};

}

/// \brief Initializes a test and checkes the command-line arguments.
///
/// If qthreads is being used, this also initializes the qthreads library.
///
/// \param argc The integer argc from main().
/// \param argv The array argv from main().
void init_test(int argc, char* argv[])
{
  bool arg_error = argc < NUM_FILE_ARGS ? true : false;

  if ((!arg_error && strcmp(argv[1], "rmat") == 0 && argc < NUM_RMAT_ARGS) ||
       (!arg_error && strcmp(argv[1], "er") == 0 && argc < NUM_ER_ARGS))
  {
    arg_error = true;
  }

  if (arg_error)
  {
#if defined(USING_QTHREADS) || defined (_OPENMP)
    fprintf(stderr, "\nSYNTAX\n"
                    "    %s { generator | filename } num_threads\n\n", argv[0]);
#else
    fprintf(stderr, "\nSYNTAX\n"
                    "    %s { generator | filename }\n\n", argv[0]);
#endif

    detail::print_test_description();

    exit(1);
  }

#if defined(USING_QTHREADS) || defined (_OPENMP)
  int nt_pos = 2;
  if (strcmp(argv[1], "rmat") == 0)
  {
    nt_pos = 3;
  }
  else if (strcmp(argv[1], "er") == 0)
  {
    nt_pos = 4;
  }

#ifdef USING_QTHREADS
  qthread_init(atoi(argv[nt_pos]));
#else
  omp_set_num_threads(atoi(argv[nt_pos]));
#endif
#endif
}

/// \brief Initializes a test and checkes the command-line arguments.
///
/// If qthreads is being used, this also initializes the qthreads library.
///
/// \param argc The integer argc from main().
/// \param argv The array argv from main().
/// \param num_ts_args The number of extra arguments specified by the test.
///                    This is the size of ts_arg_names and ts_arg_descs.
/// \param ts_arg_names The names of the extra arguments.
/// \param ts_arg_descs The descriptions of the extra arguments.
/// \param ts_argc The actual number of extra arguments entered by the user.
///                Analogous to argc in main().
/// \param ts_argv An array containing the extra arguments entered by the
///                user.  Analogous to argv in main().
void init_test(int argc, char* argv[], int num_ts_args,
               const char** ts_arg_names, const char** ts_arg_descs,
               int& ts_argc, char**& ts_argv)
{
  bool arg_error = argc < NUM_FILE_ARGS ? true : false;

  int num_internal_args = NUM_FILE_ARGS;
  if (!arg_error && strcmp(argv[1], "rmat") == 0)
  {
    if (argc < NUM_RMAT_ARGS) arg_error = true;
    num_internal_args = NUM_RMAT_ARGS;
  }
  else if (!arg_error && strcmp(argv[1], "er") == 0)
  {
    if (argc < NUM_ER_ARGS) arg_error = true;
    num_internal_args = NUM_ER_ARGS;
  }

  if (arg_error || argc < num_internal_args + num_ts_args)
  {
#if defined(USING_QTHREADS) || defined(_OPENMP)
    fprintf(stderr, "\nSYNTAX\n"
                    "    %s { generator | filename } num_threads", argv[0]);
#else
    fprintf(stderr, "\nSYNTAX\n"
                    "    %s { generator | filename }", argv[0]);
#endif

    // Print the test specific arguments with intelligent line wrapping.
    int executable_length = strlen(argv[0]);
    int line_length = executable_length + 29;
    for (int i = 0; i < num_ts_args; ++i)
    {
      int arg_length = strlen(ts_arg_names[i]);
      if (line_length + arg_length + 1 > 78)
      {
        fprintf(stderr, "\n");
        for (int i = 0; i < executable_length + 4; ++i)
        {
          fprintf(stderr, " ");
        }
        line_length = executable_length + 4;
      }

      line_length += arg_length + 1;
      fprintf(stderr, " %s", ts_arg_names[i]);
    }
    fprintf(stderr, "\n\n");

    detail::print_test_description();

    char desc_copy[MAX_DESC_LEN + 1];
    char formatted_desc[MAX_FORM_DESC_LEN];

    const char* whitespace_chars = " \f\n\r\t\v";

    fprintf(stderr, "TEST SPECIFIC ARGUMENTS\n");
    for (int i = 0; i < num_ts_args; ++i)
    {
      // Get a copy of the description.
      strncpy(desc_copy, ts_arg_descs[i], MAX_DESC_LEN);
      if (strlen(ts_arg_descs[i]) >= MAX_DESC_LEN)
      {
        desc_copy[MAX_DESC_LEN] = '\0';
      }

      // Create a nicely formatted version of the description and print it.
      strcpy(formatted_desc, "       ");
      int desc_pos = 7;
      int line_length = 7;

      char* tok = strtok(desc_copy, whitespace_chars);
      while (tok != NULL)
      {
        int tok_length = strlen(tok) + 1;
        if (line_length + tok_length > 78)
        {
          strcpy(formatted_desc + desc_pos, "\n       ");
          desc_pos += 8;
          line_length = 7;
        }

        formatted_desc[desc_pos] = ' ';
        strcpy(formatted_desc + desc_pos + 1, tok);
        desc_pos += tok_length;
        line_length += tok_length;

        tok = strtok(NULL, whitespace_chars);
      }

      fprintf(stderr, "    %s\n%s\n\n", ts_arg_names[i], formatted_desc);
    }

    exit(1);
  }

  ts_argv = argv + num_internal_args;
  ts_argc = argc - num_internal_args;

#if defined(USING_QTHREADS) || defined (_OPENMP)
  int nt_pos = 2;
  if (strcmp(argv[1], "rmat") == 0)
  {
    nt_pos = 3;
  }
  else if (strcmp(argv[1], "er") == 0)
  {
    nt_pos = 4;
  }

#ifdef USING_QTHREADS
  qthread_init(atoi(argv[nt_pos]));
#else
  omp_set_num_threads(atoi(argv[nt_pos]));
#endif
#endif
}

/// \brief Generates or reads a graph.
///
/// \param g Graph object into which the read or generated graph is placed.
/// \param argc The integer argc from main().
/// \param argv The array argv from main().
template <typename Graph>
void create_test_graph(Graph& g, int argc, char* argv[])
{
  int llen = strlen(argv[1]);

  if (strcmp(argv[1], "rmat") == 0)
  {
    generate_rmat_graph(g, atoi(argv[2]), 8);
  }
  else if (strcmp(argv[1], "er") == 0)
  {
    generate_erdos_renyi_graph(g, atoi(argv[2]), atof(argv[3]));
  }
  else if (llen > 3 && strcmp(&argv[1][llen-3], "mtx") == 0)
  {
    // Matrix-market input.
    read_matrix_market(g, argv[1]);
  }
  else if (llen > 4 && strcmp(&argv[1][llen-4], "srcs") == 0)
  {
    // Snapshot input: <file>.srcs, <file>.dests
    char* srcs_fname = argv[1];
    char dests_fname[256];

    strcpy(dests_fname, srcs_fname);
    strcpy(&dests_fname[llen - 4], "dests");

    read_binary(g, srcs_fname, dests_fname);
  }
  else if (llen > 6 && strcmp(&argv[1][llen-6], "dimacs") == 0)
  {
    // DIMACS input.
    read_dimacs(g, argv[1]);
  }
  else if (llen > 4 && strcmp(&argv[1][llen-4], "mmap") == 0)
  {
    detail::create_mmap_graph<Graph, mmap_traits<Graph>::type> cmg;
    cmg(g, argv);
  }
  else
  {
    fprintf(stderr, "Invalid filename or generate type: %s\n", argv[1]);
    exit(1);
  }

  if (num_vertices(g) == 0)
  {
    fprintf(stderr, "Error reading or creating file.  Exiting.\n");
    exit(1);
  }
}

/// \brief Generates or reads a graph including edge weights, if available.
///
/// \param g Graph object into which the read or generated graph is placed.
/// \param weights Array into which the edge weights are placed, if available.
/// \param argc The integer argc from main().
/// \param argv The array argv from main().
///
/// Note that edge weights are currently available only for R-MAT and DIMACS
/// graph formats.
template <typename Graph, typename T>
void create_test_graph(Graph& g, dynamic_array<T>& weights,
                       int argc, char* argv[])
{
  int llen = strlen(argv[1]);

  if (strcmp(argv[1], "rmat") == 0)
  {
    generate_rmat_graph(g, atoi(argv[2]), 8);
  }
  else if (strcmp(argv[1], "er") == 0)
  {
    generate_erdos_renyi_graph(g, atoi(argv[2]), atof(argv[3]));
  }
  else if (llen > 3 && strcmp(&argv[1][llen-3], "mtx") == 0)
  {
    // Matrix-market input.
    read_matrix_market(g, argv[1], weights);
  }
  else if (llen > 4 && strcmp(&argv[1][llen-4], "srcs") == 0)
  {
    // Snapshot input: <file>.srcs, <file>.dests
    char* srcs_fname = argv[1];
    char dests_fname[256];

    strcpy(dests_fname, srcs_fname);
    strcpy(&dests_fname[llen - 4], "dests");

    read_binary(g, srcs_fname, dests_fname);
  }
  else if (llen > 6 && strcmp(&argv[1][llen-6], "dimacs") == 0)
  {
    // DIMACS input.
    read_dimacs(g, argv[1], weights);
  }
  else if (llen > 4 && strcmp(&argv[1][llen-4], "mmap") == 0)
  {
    detail::create_mmap_graph<Graph, mmap_traits<Graph>::type> cmg;
    cmg(g, argv);
  }
  else
  {
    fprintf(stderr, "Invalid filename or generate type: %s\n", argv[1]);
    exit(1);
  }

  if (num_vertices(g) == 0)
  {
    fprintf(stderr, "Error reading or creating file.  Exiting.\n");
    exit(1);
  }
}

}

#undef MAX_DESC_LEN
#undef MAX_FORM_DESC_LEN
#undef NUM_RMAT_ARGS
#undef NUM_ER_ARGS

#endif
