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
/*! \file test_sssp_deltastepping.cpp

    \author William McLendon (wcmclen@sandia.gov)

    \date 1/7/2008
*/
/****************************************************************************/

#include <mtgl/util.hpp>
#include <mtgl/compressed_sparse_row_graph.hpp>
#include <mtgl/sssp_deltastepping.hpp>
#include <mtgl/mtgl_test.hpp>
#include <mtgl/random.hpp>

using namespace mtgl;

#define ASSERT 1

#ifdef __MTA__

#define BEGIN_LOGGING() \
  mta_resume_event_logging(); \
  int issues      = mta_get_task_counter(RT_ISSUES); \
  int streams     = mta_get_task_counter(RT_STREAMS); \
  int concurrency = mta_get_task_counter(RT_CONCURRENCY); \
  int traps       = mta_get_task_counter(RT_TRAP);

#define END_LOGGING() { \
  issues      = mta_get_task_counter(RT_ISSUES) - issues; \
  streams     = mta_get_task_counter(RT_STREAMS) - streams; \
  concurrency = mta_get_task_counter(RT_CONCURRENCY) - concurrency; \
  traps       = mta_get_task_counter(RT_TRAP) - traps; \
  mta_suspend_event_logging(); \
  printf("issues: %d, streams: %lf, " \
         "concurrency: %lf, " \
         "traps: %d\n", \
         issues, \
         streams / (double) issues, \
         concurrency / (double) issues, \
         traps); \
}

#else
  #define BEGIN_LOGGING() ;
  #define END_LOGGING()   ;
#endif

template <typename T>
void generate_edge_weights(T min, T max, T** realWt, int sz,
                           int max_weight = 100)
{
  T* wgts  = (T*) malloc(sizeof(T) * sz);

  double* rvals = (double*) malloc(sizeof(double) * sz);
  mt_drand48(sz, rvals);

  #pragma mta assert parallel
  for (int i = 0; i < sz; ++i) wgts[i] = rvals[i];

  free(rvals);

  *realWt = wgts;
}

template <typename T>
void generate_range(T min, T max, T** Wt, int sz)
{
  assert(max > min);

  T* wgts  = (T*) malloc(sizeof(T) * sz);

  long* rvals = (long*) malloc(sizeof(long) * sz);
  mt_lrand48(sz, rvals);

  #pragma mta assert nodep
  for (int i = 0; i < sz; ++i)
  {
    if (min == 0)
    {
      wgts[i] = rvals[i] % (max + 1);
    }
    else
    {
      wgts[i] = rvals[i] % (max - min) + min;
    }
  }

  free(rvals);

  *Wt = wgts;
}

// ***************************************************************************
//
//      DDDDD  DDDDDD DD    DDDDDD  DDDD     SSSSSS SSSSSS SSSSSS SSSSS
//      DD  DD DD     DD      DD   DD  DD    SS       SS   SS     SS  SS
//      DD  DD DDDDD  DD      DD   DDDDDD    SSSSSS   SS   SSSSS  SSSSS
//      DD  DD DD     DD      DD   DD  DD        SS   SS   SS     SS
//      DDDDD  DDDDDD DDDDDD  DD   DD  DD    SSSSSS   SS   SSSSSS SS
///
// ***************************************************************************
template <typename graph>
double* run_sssp_deltastepping(graph& ga, int vs, double* realWt,
                               double& checksum, double& time)
{
  double* result = (double*) malloc(sizeof(double) * num_vertices(ga));
  mt_timer timer;
  sssp_deltastepping<graph, int> sssp_ds(ga, vs, realWt, &checksum, result);

  timer.start();
  sssp_ds.run();
  timer.stop();
  time = timer.getElapsedSeconds();

  return result;
}

template <typename graph>
void test_sssp_deltastepping(graph& ga, int vs, int nSrcs, int nTrials = 1)
{
  double the_time = 0.0;
  double* realWt = NULL;
  int* srcs   = NULL;
  double checksum = 0.0;
  bool started = false;
  double prev_checksum = 0.0;

  generate_edge_weights<double>(0.00, 1.00, &realWt, num_edges(ga));

  mt_srand48(3324545);
  generate_range<int>(0, num_vertices(ga) - 1, &srcs, nSrcs);

  if (vs != -1)
  {
    nSrcs   = 1;
    srcs[0] = vs;
  }

  double* result = 0;

  for (int i = 0; i < nSrcs; ++i)
  {
    printf("srcs[i]: %d\n", srcs[i]);

    for (int j = 0; j < nTrials; ++j)
    {
      double time = 0.0;
      result = run_sssp_deltastepping<graph>(ga, srcs[i], realWt,
                                             checksum, time);

      the_time += time;

      if (started && checksum != prev_checksum)
      {
        fprintf(stderr, "sssp checksum error\n");
      }

      if (nTrials > 1 && j < nTrials)
      {
        free(result);
        result = 0;
      }

      started = true;
      prev_checksum = checksum;
    }

    printf("\n");

    if (nTrials == 1)
    {
      free(result);
      result = 0;
    }
  }

  printf("Avg time: %8.3lf\n", the_time / (nTrials * nSrcs));

  free(realWt);
  realWt = NULL;
  free(srcs);
  srcs   = NULL;
}

int main(int argc, char* argv[])
{
  typedef compressed_sparse_row_graph<directedS> Graph;
  typedef graph_traits<Graph>::size_type size_type;

  mt_srand48(0);

  const int num_ts_args = 1;
  const char* ts_arg_names[num_ts_args] = { "src" };
  const char* ts_arg_descs[num_ts_args] =
  {
    "If positive, the source vertex id.  If negative, the number of "
      "randomly-generated sources to try."
  };
  char** ts_argv;
  int ts_argc;

  init_test(argc, argv, num_ts_args, ts_arg_names,
            ts_arg_descs, ts_argc, ts_argv);

  Graph ga;
  create_test_graph(ga, argc, argv);

  int s_arg = atoi(ts_argv[0]);
  int vs = (s_arg >= 0) ? s_arg : -1;
  int nSrcs = (s_arg < 0) ? -s_arg : 1;

  size_type order = num_vertices(ga);
  size_type size  = num_edges(ga);

  std::cout << "degree(" << 0 << ", g) = " << out_degree(vertices(ga)[0], ga)
            << std::endl
            << "order(g) = " << order << std::endl
            << "size(g) = " << size << std::endl;

  size_type trials = 1;
  vertex_id_map<Graph> vid_map = get(_vertex_id_map, ga);

  mt_timer timer;
  timer.start();
  test_sssp_deltastepping(ga, vs, nSrcs, trials);
  timer.stop();
  printf("sssp secs: %f\n", timer.getElapsedSeconds());

  return 0;
}
