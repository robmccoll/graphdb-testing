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
/*! \file sssp_deltastepping.hpp

    \brief Deltastepping SSSP code.

    \author William McLendon (wcmclen@sandia.gov)

    \date 1/7/2008
*/
/****************************************************************************/

#ifndef MTGL_SSSP_DELTASTEPPING_HPP
#define MTGL_SSSP_DELTASTEPPING_HPP

#include <cmath>

#include <mtgl/mtgl_adapter.hpp>
#include <mtgl/visit_adj.hpp>

namespace mtgl {

#define MTGL_USE_VISITADJ   0    // enable the visit_adj() use via MTGL
#define MEASURE_TRAPS       0
#define SSSP_INIT_SIZE 100

// sort weights ==> reordering edges in the underlying graph, which is a
//                  side-effect.  FWIW, these algorithms should be side-effect
//                  free thus this option should be disabled.

// Some useful tools for debugging / printing messages
//#define PRINTDEBUG 0
//#if PRINTDEBUG
//  #define INFO(...) { printf("INFO :\t");printf(__VA_ARGS__);fflush(stdout); }
//#else
//  #define INFO(...) ;
//#endif

#ifdef DEBUG
        #define mt_assert(x) if (!x) {\
    printf("%s:%i: assertion failed: ", __FILE__, __LINE__);\
    assert(x);\
}
#endif

#define CREATE_TIMER(t)  mt_timer t;
#define START_TIMER(t)   t.start();
#define STOP_TIMER(t, tag)    {\
    t.stop();\
    printf("%25s\t%9.4f\n", tag, t.getElapsedSeconds()); fflush(stdout);\
}

#if defined(__MTA__) && MEASURE_TRAPS
        #define START_TRAP(t) t = mta_get_task_counter(RT_TRAP);
        #define STOP_TRAP(t, tag) {\
    t = mta_get_task_counter(RT_TRAP) - t;\
    printf("%25s\t%9d\n", tag, t); fflush(stdout);\
}
#else
        #define START_TRAP(t)     ;
        #define STOP_TRAP(t, tag) ;
#endif

template<typename graph_t, typename timestamp_t>
class sssp_deltastepping {
public:
  typedef typename graph_traits<graph_t>::size_type size_type;
  typedef typename graph_traits<graph_t>::vertex_descriptor vertex_t;
  typedef typename graph_traits<graph_t>::edge_descriptor edge_t;
  typedef typename graph_traits<graph_t>::out_edge_iterator EITER_T;
  typedef typename graph_traits<graph_t>::vertex_iterator VITER_T;

  // G        [ IN] - reference to graph adapter.
  // s        [ IN] - source vertex id
  // realWt   [ IN] - edge weights (|E|)
  // cs       [OUT] - pointer to checksum value (scalar)
  // result   [OUT] - vertex weights (|V|)
  // vertTime [ IN] - vertex time stamps (optional) (|V|)
  sssp_deltastepping(graph_t& _G, long _s, double* _realWt, double* _cs,
                     double* _result, timestamp_t* _vertTime = NULL) :
    G(_G), s(_s), realWt(_realWt), cs(_cs), result(_result), vertTime(_vertTime)
  {
//    printAdjacency(s);
//    printGraph();
  }

  double run(void)
  {
#if defined(__MTA__) && MEASURE_TRAPS
    long __traps;
#endif

    vertex_id_map<graph_t> vid_map = get(_vertex_id_map, G);
    edge_id_map<graph_t>   eid_map = get(_edge_id_map,   G);

    mt_timer mttimer;

#ifdef DEBUG
    mt_assert(realWt != NULL);
    mt_assert(result != NULL);
#endif

//    s = 1;
    size_type m = num_edges(G);
    size_type n = num_vertices(G);

    VITER_T verts = vertices(G);

    size_type Gn = n;

    double* Wt = realWt;

    double maxWt = 0.0;
    for (size_type i = 0; i < m; ++i)
    {
      if (Wt[i] > maxWt) maxWt = Wt[i];
    }

    double INFTY = maxWt * 20;              // 20 is a surrogate for diameter.
    double delta = ((double) n) / ((double) m);
    double delta_l = delta;
    double INFTY_l = INFTY;
    long numBuckets = (long) (INFTY / delta + 1);  //(long)(1.0/delta) * (n/2);

#ifdef DEBUG
    fprintf(stderr, "source: %ld\n", s);
    fprintf(stderr, "delta: %lf\n", delta);
    fprintf(stderr, "No. of buckets: %ld\n", numBuckets);

    // Assuming a max of n phases.  Storing the no. of verts visited in each
    // phase for statistics.
    long* numVertsInPhase = (long*) calloc(n, sizeof(long));
#endif

    mttimer.start();

    /* Memory allocation */
    #pragma mta trace "malloc start"

    long** Buckets = (long**) malloc(numBuckets * sizeof(long*));
    long*  Bsize   = (long*) calloc(numBuckets, sizeof(long));
    long*  Bcount  = (long*) calloc(numBuckets, sizeof(long));
    long*  Braise  = (long*) calloc(4 * numBuckets, sizeof(long));
    long*  Bdel    = (long*) calloc(numBuckets, sizeof(long));

    // GEM: This big chunk memory allocation needs to go.  There should not
    //      be much of a speed advantage to allocating 1 large block vs. 9
    //      smaller ones.  The headache of the allocation is not worth the
    //      miniscule time savings.  The only reason I haven't nixed it yet
    //      is that I'm concerned there is some weird overlapping of the
    //      arrays going on that isn't obvious.  Just the fact that this is
    //      something you have to be concerned with is a big red flag that
    //      this method is VERY wrong.
    #pragma mta trace "memblock malloc start"
    long* memBlock = (long*) malloc((9 * n + 1) * sizeof(long));
    #pragma mta trace "memblock malloc end"

    long* S          = memBlock;
    long* R          = memBlock + n;
    long* SLoc       = memBlock + 2 * n;
    long* RLoc       = memBlock + 3 * n;
    long* Btemp      = memBlock + 4 * n;
    long* vBucketMap = memBlock + 5 * n;
    long* vBucketLoc = memBlock + 6 * n;
    double* dR = (double*)(memBlock + 7 * n);
    double* d = (double*)(memBlock + 8 * n);

    /* Initializing the bucket data structure */
    #pragma mta assert nodep
    for (size_type i = 0; i < n; ++i)
    {
      vBucketMap[i] = -1;
      d[i]    = INFTY_l;
      RLoc[i] = -1;
      SLoc[i] = -1;
    }
    d[n] = INFTY_l;

    #pragma mta trace "malloc end"

    R[0]   = s;
    dR[0]  = 0;
    long Rcount = 0;
    long Scount = 0;

//CREATE_TIMER(__t1);
//START_TIMER(__t1);

    #pragma mta trace "relax s start"
    long lastBucketNum = relax(R, dR, RLoc, 1, G, d, Buckets, 0, numBuckets,
                               Bsize, Bcount, Braise, Bdel, vBucketMap,
                               vBucketLoc, delta_l, INFTY_l);
    #pragma mta trace "relax s end"

//STOP_TIMER(__t1, "relax s");

#ifdef DEBUG
    long numRelaxations = 1;
    long numRelaxationsHeavy = 0;
    long numPhases = 0;
    long maxBucketNum = 0;
#endif

    long currBucketNum = 0;
    while (currBucketNum <= lastBucketNum)
    {
      if (Bcount[currBucketNum] == 0)
      {
        ++currBucketNum;
        continue;
      }

      /* Handle light edges */
      while (Bcount[currBucketNum] != 0)
      {
        long Bcount_t = Bcount[currBucketNum];
        long Bdel_t = Bdel[currBucketNum];
        long* Bt = Buckets[currBucketNum];
        long Btemp_count = 0;
        Rcount = 0;

//CREATE_TIMER(__t0);
//START_TIMER(__t0);

        if (Bdel_t == Bcount_t)
        {
          Btemp_count = 0;
          /* The bucket doesn't have a lot of empty spots */
        }
        else if (Bdel_t < Bcount_t / 3 + 2)
        {
//          INFO("\tlight nodup\n");
          START_TRAP(__traps);
//CREATE_TIMER(__t2);
//START_TIMER(__t2);

          Btemp_count = Bcount_t;

          if (Bcount_t > 30)
          {
            #pragma mta trace "light nodup start"
#if MTGL_USE_VISITADJ
            #pragma mta assert parallel
            #pragma mta block dynamic schedule
            for (long i = 0; i < Bcount_t; ++i)
            {
              const size_type u = Bt[i];  // source vert
              if (u == Gn) continue;

              const double du = d[u];

              light_collect_visadj_t lnodup_vis(G, delta_l, Wt,
                                                d, du, dR, Rcount, R, RLoc,
                                                vertTime, vid_map, eid_map);

              visit_adj<graph_t, light_collect_visadj_t>(G, u, lnodup_vis, 1);
            }
#else
            #pragma mta assert parallel
//            #pragma mta block dynamic schedule
            #pragma mta loop future
            for (long i = 0; i < Bcount_t; ++i)
            {
              const size_type u = Bt[i];  // source vert
              if (u == Gn) continue;

              const double du = d[u];

              const vertex_t vu = verts[u];
              const size_type deg = out_degree(vu, G);
              EITER_T inc_edges = out_edges(vu, G);

              #pragma mta assert parallel
              for (size_type ineigh = 0; ineigh < deg; ++ineigh)
              {
                edge_t e = inc_edges[ineigh];
                const size_type j = get(eid_map, e);
                const size_type v = get(vid_map, target(e, G));

//                INFO("\tP light nodup: e[%lu]  %lu->%lu\n", j, u, v);

                if (vertTime == NULL || vertTime[u] <= vertTime[v])
                {
                  if (du + Wt[j] < d[v])
                  {
                    if (delta_l > Wt[j])
                    {
                      long rlv = mt_readfe(RLoc[v]);

                      if (rlv == -1)
                      {
                        long pos = mt_incr(Rcount, 1);
                        R[pos] = v;
                        dR[pos] = du + Wt[j];
                        mt_write(RLoc[v], pos);
                      }
                      else
                      {
                        if (du + Wt[j] < dR[rlv]) dR[rlv] = du + Wt[j];

                        mt_write(RLoc[v], rlv);
                      }
                    }
                  }
                }
              }
            }
#endif
            #pragma mta trace "light nodup end"
          }
          else
          {
#if MTGL_USE_VISITADJ
            for (long i = 0; i < Bcount_t; ++i)
            {
              size_type u = Bt[i];
              if (u == Gn) continue;

              double du = d[u];

              light_collect_visadj_t light_vis(G, delta_l, Wt,
                                               d, du, dR, Rcount, R, RLoc,
                                               vertTime, vid_map, eid_map);

              visit_adj<graph_t, light_collect_visadj_t>(G, u, light_vis, 1);
            }
#else
            for (long i = 0; i < Bcount_t; ++i)
            {
              const size_type u  = Bt[i];
              if (u == Gn) continue;

              const double du  = d[u];

              vertex_t vu  = verts[u];
              const size_type deg  = out_degree(vu, G);
              EITER_T inc_edges = out_edges(vu, G);

              #pragma mta assert parallel
              for (size_type ineigh = 0; ineigh < deg; ++ineigh)
              {
                edge_t e = inc_edges[ineigh];
                const size_type j = get(eid_map, e);
                const size_type v = get(vid_map, target(e, G));

//                INFO("\tS light nodup: e[%3lu]  %2lu->%-2lu\tj u v\n", j, u, v);
                if (vertTime == NULL || vertTime[u] <= vertTime[v])
                {
                  if (du + Wt[j] < d[v])
                  {
                    if (delta_l > Wt[j])
                    {
                      long rlv = mt_readfe(RLoc[v]);

                      if (rlv == -1)
                      {
                        long pos = mt_incr(Rcount, 1);
                        R[pos] = v;
                        dR[pos] = du + Wt[j];
                        mt_write(RLoc[v], pos);
                      }
                      else
                      {
                        if (du + Wt[j] < dR[rlv]) dR[rlv] = du + Wt[j];

                        mt_write(RLoc[v], rlv);
                      }
                    }
                  }
                }
              }
            }
#endif
          }

//STOP_TIMER(__t2, "light nodup");
          STOP_TRAP(__traps, "light nodup traps");
//CREATE_TIMER(__t3);
//START_TIMER(__t3);

          /* Add to S */
          if (Bcount_t > 30)
          {
            #pragma mta trace "light nodup S start"
            #pragma mta assert parallel
            for (long i = 0; i < Bcount_t; ++i)
            {
              size_type Gn = n;
              size_type u = Bt[i];

              if (u == Gn) continue;

              long slv = mt_readfe(SLoc[u]);

              /* Add the vertex to S */
              if (slv == -1)
              {
                long pos = mt_incr(Scount, 1);
                S[pos] = u;
                mt_write(SLoc[u], pos);
              }
              else
              {
                mt_write(SLoc[u], slv);
              }

              vBucketMap[u] = -1;
              vBucketLoc[u] = -1;
            }

            #pragma mta trace "light nodup S end"
          }
          else
          {
            for (long i = 0; i < Bcount_t; ++i)
            {
              size_type Gn = n;
              size_type u = Bt[i];

              if (u == Gn) continue;

              long slv = mt_readfe(SLoc[u]);

              /* Add the vertex to S */
              if (slv == -1)
              {
                long pos = mt_incr(Scount, 1);
                S[pos] = u;
                mt_write(SLoc[u], pos);
              }
              else
              {
                mt_write(SLoc[u], slv);
              }

              vBucketMap[u] = -1;
              vBucketLoc[u] = -1;
            }
          }
//STOP_TIMER(__t3, "light nodup S");
        }
        else
        {
          /* Bdel_t > Bcount_t/3  */
          /* There are a significant number of empty spots in the bucket.
           * So we get the non-empty vertices and store them in a compact
           * array */

          size_type Gn = n;

//CREATE_TIMER(__t4);
//START_TIMER(__t4);

          if (Bcount_t > 30)
          {
            #pragma mta trace "light dup filter start"
            #pragma mta assert nodep
            #pragma mta interleave schedule
            // #pragma mta use 60 streams
            for (long i = 0; i < Bcount_t; ++i)
            {
              size_type u = Bt[i];

              if (u != Gn)
              {
                long pos = mt_incr(Btemp_count, 1);
                Btemp[pos] = u;
              }
            }

            #pragma mta trace "light dup filter end"
          }
          else
          {
            for (long i = 0; i < Bcount_t; ++i)
            {
              size_type u = Bt[i];

              if (u != Gn)
              {
                long pos = mt_incr(Btemp_count, 1);
                Btemp[pos] = u;
              }
            }
          }
//STOP_TIMER(__t4, "light dup filter");

          /* The nested loop can be collapsed, but this results
           * in a lot of hotspots */

//CREATE_TIMER(__t5);
//START_TIMER(__t5);

//          INFO("\tlight dup\n");
          if (Btemp_count > 30)
          {
#if MTGL_USE_VISITADJ
            #pragma mta trace "light dup start"
            #pragma mta assert parallel
            for (long i = 0; i < Btemp_count; ++i)
            {
              size_type u = Btemp[i];
              double du = d[u];

              light_collect_visadj_t light_vis(G, delta_l, Wt,
                                               d, du, dR, Rcount, R, RLoc,
                                               vertTime, vid_map, eid_map);

              // visitadj in serial (note parcutoff = m = |E|).
              visit_adj<graph_t, light_collect_visadj_t>(G, u, light_vis, m);
            }
#else
            #pragma mta trace "light dup start"
            #pragma mta assert parallel
            for (long i = 0; i < Btemp_count; ++i)
            {
              const size_type u = Btemp[i];
              const double du = d[u];

              vertex_t vu = verts[u];
              const size_type deg = out_degree(vu, G);
              EITER_T inc_edges = out_edges(vu, G);

              for (size_type ineigh = 0; ineigh < deg; ++ineigh)
              {
                edge_t e = inc_edges[ineigh];
                const size_type j = get(eid_map, e);
                const size_type v = get(vid_map, target(e, G));

//                INFO("\tP light dup: e[%2lu]  %2lu->%-2lu\tj u v\n", j, u, v);

                if (vertTime == NULL || vertTime[u] <= vertTime[v])
                {
                  if (du + Wt[j] < d[v])
                  {
                    if (delta_l > Wt[j])
                    {
                      long rlv = mt_readfe(RLoc[v]);

                      if (rlv == -1)
                      {
                        long pos = mt_incr(Rcount, 1);
                        R[pos] = v;
                        dR[pos] = du + Wt[j];
                        mt_write(RLoc[v], pos);
                      }
                      else
                      {
                        if (du + Wt[j] < dR[rlv]) dR[rlv] = du + Wt[j];

                        mt_write(RLoc[v], rlv);
                      }
                    }
                  }
                }
              }
            }

            #pragma mta trace "light dup end"
#endif                                                  // MTGL_USE_VISITADJ
          }
          else
          {
#if MTGL_USE_VISITADJ
            for (long i = 0; i < Btemp_count; ++i)
            {
              size_type u = Btemp[i];
              double du = d[u];

              light_collect_visadj_t light_vis(G, delta_l, Wt,
                                               d, du, dR, Rcount, R, RLoc,
                                               vertTime, vid_map, eid_map);

              // visitadj in serial (note parcutoff = m = |E|).

              visit_adj<graph_t, light_collect_visadj_t>(G, u, light_vis, m);
            }

#else
            for (long i = 0; i < Btemp_count; ++i)
            {
              const size_type u = Btemp[i];
              const double du = d[u];

              vertex_t vu = verts[u];
              const size_type deg = out_degree(vu, G);
              EITER_T inc_edges = out_edges(vu, G);

              for (size_type ineigh = 0; ineigh < deg; ++ineigh)
              {
                edge_t e = inc_edges[ineigh];
                const size_type j = get(eid_map, e);
                const size_type v = get(vid_map, target(e, G));

//                INFO("\tP light dup: e[%2lu]  %2lu->%-2lu\tj u v\n", j, u, v);

                if (vertTime == NULL || vertTime[u] <= vertTime[v])
                {
                  if (du + Wt[j] < d[v])
                  {
                    if (delta_l > Wt[j])
                    {
                      long rlv = mt_readfe(RLoc[v]);

                      if (rlv == -1)
                      {
                        long pos = mt_incr(Rcount, 1);
                        R[pos] = v;
                        dR[pos] = du + Wt[j];
                        mt_write(RLoc[v], pos);
                      }
                      else
                      {
                        if (du + Wt[j] < dR[rlv]) dR[rlv] = du + Wt[j];

                        mt_write(RLoc[v], rlv);
                      }
                    }
                  }
                }
              }
            }
#endif                                                  // MTGL_USE_VISITADJ
          }

//STOP_TIMER(__t5, "light dup");
//CREATE_TIMER(__t6);
//START_TIMER(__t6);
//          INFO("\tlight dup S\n");

          if (Btemp_count > 30)
          {
            #pragma mta trace "light dup S start"
            #pragma mta assert parallel
            for (long i = 0; i < Btemp_count; ++i)
            {
              size_type u = Btemp[i];
              long slv = mt_readfe(SLoc[u]);

              /* Add the vertex to S */

              if (slv == -1)
              {
                long pos = mt_incr(Scount, 1);
                S[pos] = u;
                mt_write(SLoc[u], pos);
              }
              else
              {
                mt_write(SLoc[u], slv);
              }

              vBucketMap[u] = -1;
              vBucketLoc[u] = -1;
            }

            #pragma mta trace "light dup S end"
          }
          else
          {
            for (long i = 0; i < Btemp_count; ++i)
            {
              size_type u = Btemp[i];
              long slv = mt_readfe(SLoc[u]);

              /* Add the vertex to S */

              if (slv == -1)
              {
                long pos = mt_incr(Scount, 1);
                S[pos] = u;
                mt_write(SLoc[u], pos);
              }
              else
              {
                mt_write(SLoc[u], slv);
              }

              vBucketMap[u] = -1;
              vBucketLoc[u] = -1;
            }
          }

//STOP_TIMER(__t6, "light dup S");

        }        /* end of if .. then .. else */
                 /* We have collected all the light edges in R */

//STOP_TIMER(__t0, "collect edges");

#ifdef DEBUG
        if (Btemp_count != 0) { maxBucketNum = currBucketNum; }
#endif

        Bcount[currBucketNum] = 0;
        Bdel[currBucketNum] = 0;

#ifdef DEBUG
        numVertsInPhase[numPhases++] = Rcount;
        numRelaxations += Rcount;
#endif

        /* Relax all light edges */
        if (Rcount != 0)
        {
          lastBucketNum = relax(R, dR, RLoc, Rcount, G, d, Buckets,
                                currBucketNum, numBuckets, Bsize,
                                Bcount, Braise, Bdel, vBucketMap,
                                vBucketLoc, delta_l, INFTY_l);
        }
      }                         /* inner-most while end */

      Rcount = 0;

//CREATE_TIMER(__hvy);
//START_TIMER(__hvy);

      /* Collect heavy edges into R */
//      INFO("\theavy collect\n");

      if (Scount > 10)
      {
#if MTGL_USE_VISITADJ
        #pragma mta trace "heavy collect start"
        #pragma mta assert parallel
        #pragma mta block dynamic schedule
        for (long i = 0; i < Scount; ++i)
        {
          size_type u = S[i];
          SLoc[u] = -1;
          double du = d[u];

          heavy_collect_visadj_t hc_vis(G, delta_l, Wt, d, du, dR,
                                        Rcount, R, RLoc, vertTime,
                                        vid_map, eid_map);

          visit_adj<graph_t, heavy_collect_visadj_t>(G, u, hc_vis, m);
        }
#else

        #pragma mta trace "heavy collect start"
        #pragma mta assert parallel
        #pragma mta block dynamic schedule
        for (long i = 0; i < Scount; ++i)
        {
          const size_type u = S[i];
          SLoc[u] = -1;
          const double du = d[u];

          vertex_t vu = verts[u];
          const size_type deg = out_degree(vu, G);
          EITER_T inc_edges = out_edges(vu, G);

          for (size_type ineigh = 0; ineigh < deg; ++ineigh)
          {
            edge_t e = inc_edges[ineigh];
            const size_type j = get(eid_map, e);
            const size_type v = get(vid_map, target(e, G));

//            INFO("\tP heavy collect: e[%2lu] %2lu->%-2lu  j u v\n", j, u, v);

            if (vertTime == NULL || vertTime[u] <= vertTime[v])
            {
              if (Wt[j] + du < d[v])
              {
                if (delta_l <= Wt[j])  // SORTW
                {
                  const long rlv = mt_readfe(RLoc[v]);

                  if (rlv == -1)
                  {
                    long pos = mt_incr(Rcount, 1);
                    R[pos]  = v;
                    dR[pos] = d[u] + Wt[j];
                    mt_write(RLoc[v], pos);
                  }
                  else
                  {
                    if (du + Wt[j] < dR[rlv]) dR[rlv] = du + Wt[j];

                    mt_write(RLoc[v], rlv);
                  }
                }
              }
            }
          }
        }

        #pragma mta trace "heavy collect end"
#endif
      }
      else
      {
#if MTGL_USE_VISITADJ
        for (long i = 0; i < Scount; ++i)
        {
          size_type u = S[i];
          SLoc[u] = -1;
          double du = d[u];

          heavy_collect_visadj_t hc_vis(G, delta_l, Wt, d, du, dR,
                                        Rcount, R, RLoc, vertTime,
                                        vid_map, eid_map);

          visit_adj<graph_t, heavy_collect_visadj_t>(G, u, hc_vis, m);
        }
#else
        for (long i = 0; i < Scount; ++i)
        {
          const size_type u = S[i];
          SLoc[u] = -1;
          const double du = d[u];

          vertex_t vu = verts[u];
          const size_type deg = out_degree(vu, G);
          EITER_T inc_edges = out_edges(vu, G);

          for (size_type ineigh = 0; ineigh < deg; ++ineigh)
          {
            edge_t e = inc_edges[ineigh];
            const size_type j = get(eid_map, e);
            const size_type v = get(vid_map, target(e, G));

//            INFO("\tS heavy collect: e[%2lu] %2lu->%-2lu  j u v\n", j, u, v);

            if (vertTime == NULL || vertTime[u] <= vertTime[v])
            {
              if (Wt[j] + du < d[v])
              {
                if (delta_l <= Wt[j])  // SORTW
                {
                  const long rlv = mt_readfe(RLoc[v]);

                  if (rlv == -1)
                  {
                    const long pos = mt_incr(Rcount, 1);
                    R[pos]  = v;
                    dR[pos] = d[u] + Wt[j];
                    mt_write(RLoc[v], pos);
                  }
                  else
                  {
                    if (du + Wt[j] < dR[rlv]) dR[rlv] = du + Wt[j];

                    mt_write(RLoc[v], rlv);
                  }
                }
              }
            }
          }
        }
#endif
      }
//STOP_TIMER(__hvy, "heavy collect");
      Scount = 0;

      /* Relax heavy edges */

      #pragma mta trace "heavy relax start"
      if (Rcount != 0)
      {
//CREATE_TIMER(__tmr);
//START_TIMER(__tmr);

        lastBucketNum = relax(R, dR, RLoc, Rcount, G, d, Buckets,
                              currBucketNum, numBuckets, Bsize, Bcount,
                              Braise, Bdel, vBucketMap, vBucketLoc,
                              delta_l, INFTY_l);

//STOP_TIMER(__tmr, "heavy relax");
      }

      #pragma mta trace "heavy relax end"

#ifdef DEBUG
      numRelaxationsHeavy += Rcount;
#endif
    }

    mttimer.stop();
    double running_time = mttimer.getElapsedSeconds();

    /* Compute d checksum  AND save vert weights into result */
    double checksum = 0;
    size_type not_connected = 0;

    #pragma mta assert parallel
    for (size_type i = 0; i < n; ++i)
    {
      result[i] = d[i];   // Save Result.

      if (d[i] < INFTY_l)
      {
        checksum = checksum + d[i];
      }
      else
      {
        mt_incr(not_connected, 1);
      }
    }

    *cs = checksum;

#ifdef DEBUG
    /* Compute W checksum */
    double Wsum = 0;
    for (size_type i = 0; i < m; ++i) Wsum = Wsum + Wt[i];

    fprintf(stderr, "d checksum: %lf, W checksum: %lf, Avg. distance %lf\n",
            checksum, Wsum, checksum / (n - not_connected));
    fprintf(stderr, "Last non-empty bucket: %d\n", maxBucketNum);
    fprintf(stderr, "No. of phases: %d\n", numPhases);
    fprintf(stderr, "No. of light relaxations: %d\n", numRelaxations);
    fprintf(stderr, "No. of heavy relaxations: %d\n", numRelaxationsHeavy);
    fprintf(stderr, "Avg. no. of light edges relaxed in a phase: %d\n",
            numRelaxations / numPhases);

    if (maxBucketNum != 0)
    {
      fprintf(stderr, "Avg. no. of heavy edges relaxed per bucket: %d\n",
              numRelaxationsHeavy / maxBucketNum);
    }

    fprintf(stderr, "Total no. of relaxations: %d\n\n",
            numRelaxations + numRelaxationsHeavy);

    fflush(stderr);
#endif

    /* Free memory */
    #pragma mta assert parallel
    for (long i = 0; i < numBuckets; ++i)
    {
      if (Bsize[i] != 0) free(Buckets[i]);
    }

    free(Buckets);
    free(Bsize);
    free(Bcount);
    free(Braise);
    free(Bdel);
    free(memBlock);
#ifdef DEBUG
    free(numVertsInPhase);
#endif

//    INFO("\tChecksum = %6.1lf\n", checksum);
//    INFO("<<<\tsssp_deltastepping()\n");

    return running_time;
  }

  //-----[ end run() ]--------------------------------------------------

  //-----[ relax() ]----------------------------------------------------
  #pragma mta inline
  long relax(long* R,
             double* dR,
             long* RLoc,
             long Rcount,
             graph_t& G,
             double* d,
             long** Buckets,
             long currBucketNum,
             long numBuckets,
             long* Bsize,
             long* Bcount,
             long* Braise,
             long* Bdel,
             long* vBucketMap,
             long* vBucketLoc,
             double delta,
             double INFTY)
  {
    double delta_l = delta;
    double INFTY_l = INFTY;
    size_type Gn = num_vertices(G);
    long lastBucketNum = -1;

//    INFO(">>>\trelax()\n");

    if (Rcount > 30)
    {
      #pragma mta trace "loop 1 start"
      #pragma mta assert nodep
      #pragma mta interleave schedule
      for (long i = 0; i < Rcount; ++i)
      {
        long v = R[i];
        long bn = (long) floor(dR[i] / delta_l);
        long bn_old = vBucketMap[v];
        long offset = numBuckets * (i & 3);

        if (bn > lastBucketNum) lastBucketNum = bn;
        if (bn >= numBuckets) bn = numBuckets - 1;

/*
        if (bn >= numBuckets)
        {
          fprintf(stderr, "Error: relaxation failed, bn: %d, numBuckets: %d\n",
                  bn, numBuckets);
          exit(1);
        }
*/

        RLoc[v] = (i & 3) * Gn + mt_incr(Braise[bn + offset], 1);

        if ((d[v] < INFTY_l) && (bn_old != -1))
        {
          Buckets[bn_old][vBucketLoc[v]] = Gn;
          Bdel[bn_old]++;
        }
      }

      #pragma mta trace "loop 1 end"
    }
    else
    {
      for (long i = 0; i < Rcount; ++i)
      {
        long v = R[i];
        long bn = (long) floor(dR[i] / delta_l);
        long bn_old = vBucketMap[v];
        long offset = numBuckets * (i & 3);

        if (bn > lastBucketNum) lastBucketNum = bn;
        if (bn >= numBuckets) bn = numBuckets - 1;

        RLoc[v] = (i & 3) * Gn + mt_incr(Braise[bn + offset], 1);

        if ((d[v] < INFTY_l) && (bn_old != -1))
        {
          Buckets[bn_old][vBucketLoc[v]] = Gn;
          ++Bdel[bn_old];
        }
      }
    }

    ++lastBucketNum;

//     fprintf(stderr, "[%d %d] ", currBucketNum, lastBucketNum);

    #pragma mta trace "loop 2 start"
    for (long i = currBucketNum; i < lastBucketNum; ++i)
    {
      long* Bi = Buckets[i];
      long Bsize_i = Bsize[i];
      long size_incr = Braise[i] + Braise[i + 1 * numBuckets] +
                       Braise[i + 2 * numBuckets] + Braise[i + 3 * numBuckets];

      if ((size_incr > 0) && (Bcount[i] + size_incr >= Bsize[i]))
      {
        long Bsize_tmp = Bcount[i] + size_incr + SSSP_INIT_SIZE;
        long* Bt = (long*) malloc(Bsize_tmp * sizeof(long));

        if (Bsize_i != 0)
        {
          if (Bsize_i > 30)
          {
            #pragma mta assert nodep
            #pragma mta interleave schedule
            // #pragma mta use 60 streams
            for (long j = 0; j < Bsize_i; ++j) Bt[j] = Bi[j];
          }
          else
          {
            for (long j = 0; j < Bsize_i; ++j) Bt[j] = Bi[j];
          }

          free(Bi);
        }

        Buckets[i] = Bt;
        Bsize[i] = Bsize_tmp;
      }
    }

    #pragma mta trace "loop 2 end"

    if (Rcount > 30)
    {
      #pragma mta trace "loop 3 start"
      #pragma mta assert nodep
      #pragma mta interleave schedule
      for (long i = 0; i < Rcount; ++i)
      {
        long v = R[i];
        double x = dR[i];
        long loc = RLoc[v];
        long locDiv = loc / Gn;
        long locMod = loc % Gn;
        long bn = (long) floor(x / delta_l);
        long pos = Bcount[bn] + locMod;

        pos += (locDiv >= 1) * Braise[bn];
        pos += (locDiv >= 2) * Braise[bn + numBuckets];
        pos += (locDiv >= 3) * Braise[bn + 2 * numBuckets];

        Buckets[bn][pos] = v;
        vBucketLoc[v] = pos;
        vBucketMap[v] = bn;

//        if (d[v] != x)
//        {
//          printf("update D[%5d] <= %2.8f  (%-2.8f)\n", v, x, d[v]);
//        }

        d[v] = x;
        RLoc[v] = -1;
      }

      #pragma mta trace "loop 3 end"
    }
    else
    {
      for (long i = 0; i < Rcount; ++i)
      {
        long v = R[i];
        double x = dR[i];
        long loc = RLoc[v];
        long locDiv = loc / Gn;
        long locMod = loc % Gn;
        long bn = (long) floor(x / delta_l);
        long pos = Bcount[bn] + locMod;

        pos += (locDiv >= 1) * Braise[bn];
        pos += (locDiv >= 2) * Braise[bn + numBuckets];
        pos += (locDiv >= 3) * Braise[bn + 2 * numBuckets];

        Buckets[bn][pos] = v;
        vBucketLoc[v] = pos;
        vBucketMap[v] = bn;

//        if (d[v] != x)
//        {
//          printf("update D[%5d] <= %2.8f  (%-2.8f)\n", v, x, d[v]);
//        }

        d[v] = x;
        RLoc[v] = -1;
      }
    }

    if (lastBucketNum - currBucketNum > 30)
    {
      #pragma mta trace "loop 4 start"
      #pragma mta parallel single processor
      #pragma mta assert nodep
      #pragma mta interleave schedule
      for (long i = currBucketNum; i < lastBucketNum; ++i)
      {
        Bcount[i] += Braise[i] + Braise[i + numBuckets] +
                     Braise[i + 2 * numBuckets] + Braise[i + 3 * numBuckets];
        Braise[i] = 0;
        Braise[i + numBuckets] = 0;
        Braise[i + 2 * numBuckets] = 0;
        Braise[i + 3 * numBuckets] = 0;
      }

      #pragma mta trace "loop 4 end"
    }
    else
    {
      for (long i = currBucketNum; i < lastBucketNum; ++i)
      {
        Bcount[i] += Braise[i] + Braise[i + numBuckets] +
                     Braise[i + 2 * numBuckets] + Braise[i + 3 * numBuckets];
        Braise[i] = 0;
        Braise[i + numBuckets]   = 0;
        Braise[i + 2 * numBuckets] = 0;
        Braise[i + 3 * numBuckets] = 0;
      }
    }

//    INFO("<<<\trelax()\n");
    return lastBucketNum;
  }

  //-----[ end relax() ]------------------------------------------------

protected:
  void printGraph()
  {
    for (size_type v = 0; v < num_vertices(G); ++v) printAdjacency(v);
  }

  void printAdjacency(size_type vid)
  {
    vertex_id_map<graph_t> vid_map = get(_vertex_id_map, G);

    vertex_t vu = vertices(G)[vid];
    size_type deg = out_degree(vu, G);
    EITER_T inc_edges = out_edges(vu, G);

    printf("%lu [%lu]\t{ ", vid, deg); fflush(stdout);

    for (size_type ineigh = 0; ineigh < deg; ++ineigh)
    {
      edge_t e = inc_edges[ineigh];
      size_type tgt = get(vid_map, target(e, G));

      printf("%lu ", tgt); fflush(stdout);
    }

    printf("}\n"); fflush(stdout);
  }

  void print_double_array(double* A, size_type N, char* TAG)
  {
    for (size_type i = 0, c = 1; i < N; ++i, ++c)
    {
      printf("%6.3f ", A[i]);

      if (c == 10)
      {
        printf("\n");
        c = 0;
      }
    }

    printf("\n");
  }

  //-----[ heavy_collect_visadj_t ]-------------------------------------
  class heavy_collect_visadj_t {
  public:
    heavy_collect_visadj_t(graph_t& __ga,
                           double __delta_l,
                           double* __Wt,
                           double* __d,
                           double __du,
                           double* __dR,
                           long& __Rcount,
                           long* __R,
                           long* __RLoc,
                           timestamp_t* __vertTime,
                           vertex_id_map<graph_t>& __vid_map,
                           edge_id_map<graph_t>& __eid_map) :
      ga(__ga), delta_l(__delta_l), Wt(__Wt), d(__d),
      du(__du), dR(__dR), Rcount(__Rcount), R(__R), RLoc(__RLoc),
      vertTime(__vertTime), vid_map(__vid_map), eid_map(__eid_map) {}

    bool visit_test(vertex_t v) { return(true); }

    void operator()(edge_t& e, vertex_t& src, vertex_t& dest)
    {
      vertex_t first = source(e, ga);

      if (src == first)
      {
        size_type u = get(vid_map, src);      // source-ID
        size_type v = get(vid_map, dest);     // dest-ID
        size_type j = get(eid_map, e);        // edge-ID

//        INFO("\t- heavy collect: e[%2lu] %2lu->%-2lu   j u v\n", j, u, v);

        if (vertTime == NULL || vertTime[u] <= vertTime[v])
        {
          if (Wt[j] + du < d[v])
          {
            if (delta_l <= Wt[j])  // SORTW
            {
              long rlv = mt_readfe(RLoc[v]);

              if (rlv == -1)
              {
                long pos = mt_incr(Rcount, 1);
                R[pos]  = v;
                dR[pos] = d[u] + Wt[j];
                mt_write(RLoc[v], pos);
              }
              else
              {
                if (du + Wt[j] < dR[rlv]) dR[rlv] = du + Wt[j];

                mt_write(RLoc[v], rlv);
              }
            }
          }
        }
      }
    }

  private:
    graph_t& ga;
    timestamp_t* vertTime;
    double delta_l;
    double* Wt;
    double* d;
    double du;
    double* dR;
    long& Rcount;
    long* R;
    long* RLoc;
    vertex_id_map<graph_t>& vid_map;
    edge_id_map<graph_t>& eid_map;
  };
  //-----[ heavy_collect_visadj_t ]-------------------------------------

  //-----[ light_collect_visadj_t ]-------------------------------------
  class light_collect_visadj_t {
  public:
    light_collect_visadj_t(graph_t& __ga,
                           double __delta_l,
                           double* __Wt,
                           double* __d,
                           double __du,
                           double* __dR,
                           long& __Rcount,
                           long* __R,
                           long* __RLoc,
                           timestamp_t* __vertTime,
                           vertex_id_map<graph_t>& __vid_map,
                           edge_id_map<graph_t>& __eid_map) :
      ga(__ga), delta_l(__delta_l), Wt(__Wt), d(__d), du(__du),
      dR(__dR), Rcount(__Rcount), R(__R), RLoc(__RLoc),
      vertTime(__vertTime), vid_map(__vid_map), eid_map(__eid_map) {}

    bool visit_test(vertex_t v) { return(true); }

    void operator()(edge_t& e, vertex_t& src, vertex_t& dest)
    {
      vertex_t first = source(e, ga);

      if (src == first)
      {
        size_type u = get(vid_map, src);    // source-ID
        size_type v = get(vid_map, dest);   // dest-ID
        size_type j = get(eid_map, e);      // edge-ID

//        INFO("\tlight collect: e[%lu]  %lu->%lu\n", j, u, v);

        if (vertTime == NULL || vertTime[u] <= vertTime[v])
        {
          if (du + Wt[j] < d[v])
          {
            if (delta_l > Wt[j])  // SORTW
            {
              long rlv = mt_readfe(RLoc[v]);

              if (rlv == -1)
              {
                long pos = mt_incr(Rcount, 1);
                R[pos]  = v;
                dR[pos] = du + Wt[j];
                mt_write(RLoc[v], pos);
              }
              else
              {
                if (du + Wt[j] < dR[rlv]) dR[rlv] = du + Wt[j];

                mt_write(RLoc[v], rlv);
              }
            }
          }
        }
      }
    }

  private:
    graph_t& ga;
    timestamp_t* vertTime;
    double delta_l;
    double* Wt;
    double* d;
    double du;
    double* dR;
    long& Rcount;
    long* R;
    long* RLoc;
    vertex_id_map<graph_t>& vid_map;
    edge_id_map<graph_t>& eid_map;
  };
  //-----[ light_collect_visadj_t ]-------------------------------------

private:
  graph_t& G;
  long s;
  double* realWt;
  double* cs;
  double* result;
  timestamp_t* vertTime;
};

// clean up #defined values that should be local to this header file.

#undef MTGL_USE_VISITADJ
#undef MEASURE_TRAPS
#undef SSSP_INIT_SIZE

//#undef INFO
//#undef PRINTDEBUG

}

#endif
