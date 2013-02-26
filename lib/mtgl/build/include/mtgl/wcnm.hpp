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

/////////////////////////////////////////////////////////////////////////
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
/////////////////////////////////////////////////////////////////////////
// Author        : Aaron Clauset  (aaron@cs.unm.edu)                   //
// Location      : U. Michigan, U. New Mexico                          //
// Time          : January-August 2004                                 //
// Collaborators : Dr. Cris Moore (moore@cs.unm.edu)                   //
//                 Dr. Mark Newman (mejn@umich.edu)                    //
/////////////////////////////////////////////////////////////////////////
// --- DEEPER DESCRIPTION ---------------------------------------------
// see http://www.arxiv.org/abs/cond-mat/0408187 for more information
//

/****************************************************************************/
/*! \file wcnm.hpp

    \brief This is a weighted implementation of the CNM community detection
           algorithm.

    \author Jon Berry (jberry@sandia.gov)

    This is a slightly modified version of the FastCommunityMH code by
    Aaron Clauset.

    Note: This algorithm currently runs only in serial.  We intend to
          parallelize it later.
*/
/****************************************************************************/

#ifndef MTGL_WCNM_HPP
#define MTGL_WCNM_HPP

#define MB_THRESH 1.0

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <list>

#include <mtgl/maxheap.h>
#include <mtgl/vektor.h>

#include <mtgl/mtgl_adapter.hpp>
#include <mtgl/xmt_hash_table.hpp>
#include <mtgl/neighborhoods.hpp>

#define WRITE_JOINS_FILE 0
//#define WRITE_WEIGHTED_DIMACS 1

namespace mtgl {

void gen_merges(int nodes, char* fname, int* id, int* leader, int& npairs);
void gen_merges(int nodes, int* id, int* leader, int& npairs);
void read_merges(int nodes, char* fname, int* id, int* leader, int& npairs);

// ---------------------------------------------------------------------------
/// Edge object - defined by a pair of vertex indices and *edge pointer to next
///               in linked-list
template <typename size_type>
class Edge {
public:
  int so;              ///< Originating node.
  int si;              ///< Terminating node.
  size_type eid;       ///< MTGL edge id.
  Edge* next;          ///< Pointer for linked list of edges.

  Edge() : so(0), si(0), next(NULL) {}
  ~Edge() {}
};

// ---------------------------------------------------------------------------
/// Nodenub object - defined by a *node pointer and *node pointer.
struct nodenub {
  tuple_o* heap_ptr;   ///< Pointer to node(max,i,j) in max-heap of row maxes.
  vektor* v;           ///< Pointer stored vector for (i,j).
};

// ---------------------------------------------------------------------------
/// Tuple object - defined by an real value and (row,col) indices.
#if !defined(TUPLE_INCLUDED)
#define TUPLE_INCLUDED
struct tuple_o {
  double m;            ///< Stored value.
  double sv_diff;      ///< Support variance incr.
  int i;               ///< Row index.
  int j;               ///< Column index.
  int k;               ///< Heap index.
};
#endif

/// Ordered pair structures (handy in the program).
struct apair {
  int x;
  int y;
};

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

// ---------------------------------------------------------------------------
/// List object - simple linked list of integers.
class cnmlist {
public:
  int index;           ///< Node index.
  cnmlist* next;       ///< Pointer to next element in linked list.

  cnmlist() : index(0), next(NULL) {}
  ~cnmlist() {}
};

// ---------------------------------------------------------------------------
/// Community stub object - stub for a community list.
class stub {
public:
  bool one;
  bool valid;          ///< Is this community valid?
  int size;            ///< Size of community.
  cnmlist* members;    ///< Pointer to list of community members.
  cnmlist* last;       ///< Pointer to end of list.

  stub() : valid(false), size(0), members(NULL), last(NULL) {}

  ~stub()
  {
    cnmlist* current;
    if (members != NULL)
    {
      current = members;
      while (current != NULL)
      {
        members = current->next;
        delete current;
        current = members;
      }
    }
  }
};

// ---------------------------------------------------------------------------
// FUNCTION DECLARATIONS -----------------------------------------------------

template <typename etyp>
void buildDeltaQMatrix(bool cutstep, xmt_hash_table<int64_t, double>& ewgt,
                       etyp unit);
// JWB

void buildFilenames();
void dqSupport();
// void groupListsSetup();
void groupListsStats();
void groupListsUpdate(const int x, const int y);
void mergeCommunities(int i, int j);
bool parseCommandLine(int argc, char* argv[]);
// void readInputFile(Graph& ga);
void recordGroupLists();
void recordNetwork();

// ---------------------------------------------------------------------------
// PROGRAM PARAMETERS --------------------------------------------------------

struct netparameters {
  int n;               ///< Number of nodes in network.
  int m;               ///< Number of edges in network.
  int maxid;           ///< Maximum node id.
  int minid;           ///< Minimum node id.
  double T;            ///< Sum of edge weights.
};
netparameters gparm;

struct groupstats {
  int numgroups;       ///< Number of groups.
  double meansize;     ///< Mean size of groups.
  int maxsize;         ///< Size of largest group.
  int minsize;         ///< Size of smallest group.
  double* sizehist;    ///< Distribution of sizes.
};
groupstats gstats;

struct outparameters {
  short int textFlag;   ///< 0: no console output.
                        ///< 1: writes file outputs.
  bool suppFlag;        ///< T: no support(t) file.
                        ///< F: yes support(t) file.
  short int fileFlag;
  std::string filename; ///< Name of input file.
  std::string d_in;     ///< (dir ) directory for input file.
  std::string d_out;    ///< (dir ) directory for output file.
  std::string f_parm;   ///< (file) parameters output.
  std::string f_input;  ///< (file) input data file.
  std::string f_joins;  ///< (file) community hierarchy.
  std::string f_support;///< (file) dQ support as a function of time.
  std::string f_net;    ///< (file) .wpairs file for .cutstep network.
  std::string f_group;  ///< (file) .list of indices in communities at .cutstep.
  std::string f_gstats; ///< (file) distribution of community sizes at .cutstep.
  std::string s_label;  ///< (temp) text label for run.
  std::string s_scratch;///< (temp) text for building filenames.
  int timer;            ///< Timer for displaying progress reports.
  bool timerFlag;       ///< Flag for setting timer.
  int cutstep;          ///< Step at which to record aglomerated network.
};
outparameters ioparm;

// ---------------------------------------------------------------------------
// GLOBAL VARIABLES ----------------------------------------------------------

char pauseme;
nodenub* dq;           ///< dQ matrix.
maxheap* h;            ///< Heap of values from max_i{dQ_ij}.
double* Q;             ///< Q(t).
double* SV;            ///< Support_variance(t).
dpair Qmax;            ///< Maximum Q value and the corresponding time t.
dpair SVmin;           ///< Maximum SV value and the corresponding time t.
double* a;             ///< A_i.
double* num_nodes;     ///< #nodes per community (for WT in maxheap)
                       ///< (or a values (hence double)).
apair* joins;          ///< List of joins.
stub* c;               ///< Link-lists for communities.
double W;              ///< Sum of all edge weights.

enum {NONE};

int supportTot;
double supportAve;

// ---------------------------------------------------------------------------
// FUNCTION DEFINITIONS -------------------------------------------------------

template <typename Graph>
void readInputFile(Graph& ga,
                   Edge<typename graph_traits<Graph>::size_type>*& e,
                   Edge<typename graph_traits<Graph>::size_type>*& elist,
                   double* e_wgt, typename graph_traits<Graph>::size_type unit)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::out_edge_iterator out_edge_iterator;
  typedef typename graph_traits<Graph>::adjacency_iterator adjacency_iterator;

  vertex_id_map<Graph> vipm = get(_vertex_id_map, ga);
  edge_id_map<Graph> eipm = get(_edge_id_map, ga);

  int numnodes = 0;
  int numlinks = 0;

  time_t t1;
  t1 = time(&t1);
  time_t t2;
  t2 = time(&t2);

  size_type order = num_vertices(ga);

#ifdef DEBUG
  std::cout << " scanning input file for edge and node counts." << std::endl;
  std::cout << "  edgecount: [0]" << std::endl;
#endif

  // First scan through the input file to discover the largest node id. We
  // need to do this so that we can allocate a properly sized array for the
  // sparse matrix representation.
  vertex_iterator verts = vertices(ga);
  for (size_type i = 0; i < order; ++i)
  {
    vertex_descriptor u = verts[i];
    int vid = get(vipm, u);
    adjacency_iterator adjs = adjacent_vertices(u, ga);
    size_type deg = degree(u, ga);
    for (size_type j = 0; j < deg; ++j)
    {
      vertex_descriptor v = adjs[j];

      int s = get(vipm, u);
      int f = get(vipm, v);
      if (vid == s && (s < f))
      {
        // Count number of edges.
        ++numlinks;

        // Guarantee s < f.
        if (f < s)
        {
          int t = s;
          s = f;
          f = t;
        }

        // Track largest node index.
        if (f > numnodes) numnodes = f;

        // Check timer; if necessarsy, display.
        if (t2 - t1 > ioparm.timer)
        {
#ifdef DEBUG
          std::cout << "  edgecount: [" << numlinks << "]" << std::endl;
#endif
          t1 = t2;
          ioparm.timerFlag = true;
        }

        t2 = time(&t2);
      }
    }
  }

#ifdef DEBUG
  std::cout << "  edgecount: [" << numlinks << "] total (first pass)"
            << std::endl;
#endif

  // Store maximum index.
  gparm.maxid = numnodes + 2;

  // Create requisite number of edges.
  elist = new Edge<size_type> [2 * numlinks];

  // Index of next edge of elist to be used.
  int ecounter = 0;

  // Now that we know numnodes, we can allocate the space for the sparse
  // matrix, and then reparse the file, adding edges as necessary.
#ifdef DEBUG
  std::cout << " allocating space for network." << std::endl;
#endif

  // (unordered) sparse adjacency matrix.
  e = new Edge<size_type> [gparm.maxid + 1];

  // List of pointers to the last edge in each row.
  Edge<size_type>** last = new Edge<size_type>* [gparm.maxid + 1];

  // Numnodes now counts number of actual used node ids.
  numnodes = 0;

  // Numlinks now counts number of bi-directional edges created.
  numlinks = 0;

  // Reset timer.
  ioparm.timerFlag = false;

#ifdef DEBUG
  std::cout << " reparsing the input file to build network data structure."
            << std::endl << "  edgecount: [0]" << std::endl;
#endif

  // Sum of edge weights (from Newman's Q^w).
  double T = 0.0;

  Edge<size_type>* newedge;
  Edge<size_type>* current;     // Pointer for checking edge existence.
  bool existsFlag;              // Flag for edge existence.

  verts = vertices(ga);
  for (size_type i = 0; i < order; ++i)
  {
    vertex_descriptor u = verts[i];
    int vid = get(vipm, u);

    adjacency_iterator adjs = adjacent_vertices(u, ga);
    out_edge_iterator edgs = out_edges(u, ga);
    size_type deg = out_degree(u, ga);

    for (size_type j = 0; j < deg; ++j)
    {
      edge_descriptor edg = edgs[j];
      size_type eid = get(eipm, edg);
      vertex_descriptor v = target(edg, ga);
      int s = get(vipm, u);
      int f = get(vipm, v);

      if (vid == s && (s < f))
      {
        // Increment s,f to prevent using e[0].
        ++s;
        ++f;

        // Guarantee s < f.
        if (f < s)
        {
          int t = s;
          s = f;
          f = t;
        }

        double wgt = e_wgt[eid];
        T += wgt;

        // Increment link count (preemptive).
        ++numlinks;

        // If first edge with s, add s and (s,f).
        if (e[s].so == 0)
        {
          e[s].so = s;
          e[s].si = f;
          e[s].eid = eid;

          // Point last[s] at self.
          last[s] = &e[s];

          ++numnodes;
        }
        else                          // Try to add (s,f) to s-edgelist.
        {
          current = &e[s];
          existsFlag = false;

          // Check if (s,f) already in edgelist.
          while (current != NULL)
          {
            if (current->si == f)
            {
              existsFlag = true;      // Link already exists.
              --numlinks;             // Adjust link-count downward.
              break;
            }

            current = current->next;  // Look at next edge.
          }

          // If not already exists, append it.
          if (!existsFlag)
          {
            newedge = &elist[ecounter++];    // Grab next-free-edge.
            newedge->so = s;
            newedge->si = f;
            newedge->eid = eid;
            last[s]->next = newedge;         // Append newedge to [s]'s list.
            last[s] = newedge;               // Point last[s] to newedge.
          }
        }

        // If first edge with f, add f and (f,s).
        if (e[f].so == 0)
        {
          e[f].so = f;
          e[f].si = s;
          e[f].eid = eid;
          last[f] = &e[f];     // Point last[s] at self.
          ++numnodes;
        }
        else                               // Try to add (f,s) to f-edgelist.
        {
          // If (s,f) wasn't in s-edgelist, then.
          if (!existsFlag)
          {
            newedge = &elist[ecounter++];  // (f,s) not in f-edgelist.
            newedge->so = f;
            newedge->si = s;
            newedge->eid = eid;
            last[f]->next = newedge;       // Append newedge to [f]'s list.
            last[f] = newedge;             // Point last[f] to newedge.
          }
        }

        // Reset existsFlag.
        existsFlag = false;

        // Check timer; if necessarsy, display.
        if (t2 - t1 > ioparm.timer)
        {
#ifdef DEBUG
          std::cout << "  edgecount: [" << numlinks << "]" << std::endl;
#endif
          t1 = t2;
          ioparm.timerFlag = true;
        }

        t2 = time(&t2);
      }
    }
  }

#ifdef DEBUG
  std::cout << "  edgecount: [" << numlinks << "] total (second pass)"
            << std::endl;
#endif

  // Now, we record our work in the parameters file, and exit.
  std::ofstream fout(ioparm.f_parm.c_str(), std::ios::app);
  fout << "---NET_STATS----\n";
  fout << "MAXID-----:\t" << gparm.maxid - 2 << "\n";
  fout << "NUMNODES--:\t" << numnodes << "\n";
  fout << "NUMEDGES--:\t" << numlinks << "\n";
  fout.close();

  gparm.T = T;          // Store sum of edge weights.
  gparm.m = numlinks;   // Store actual number of edges created.
  gparm.n = numnodes;   // Store actual number of nodes used.

  delete [] last;

  return;
}

template <typename Graph>
double fastcommunity_mh(Graph& ga, double* e_wgt,
                        typename graph_traits<Graph>::size_type* leader,
                        typename graph_traits<Graph>::size_type unit,
                        int& num_communities, char* lf = 0, bool cutstep = 0)
{
  typedef typename graph_traits<Graph>::size_type size_type;

  // Default values for parameters which may be modified from the commandline.
  ioparm.timer = 20;
  ioparm.fileFlag = NONE;
  ioparm.suppFlag = false;
  ioparm.textFlag = 0;
  ioparm.s_scratch = "community";
  ioparm.s_label = "a";

  time_t t1;
  t1 = time(&t1);
  time_t t2;
  t2 = time(&t2);

  Edge<size_type>* e;      // Initial adjacency matrix (sparse).
  Edge<size_type>* elist;  // List of edges for building adjacency matrix.

  size_type order = num_vertices(ga);

  // ----------------------------------------------------------------------
  // Parse the command line, build filenames and, import the .pairs file.

#ifdef DEBUG
  std::cout << "\nFast Community Inference.\n";
  std::cout << "Copyright (c) 2004 by Aaron Clauset (aaron@cs.unm.edu)\n";
#endif

/*
   if (!parseCommandLine(argc, argv)) return 0;

   // Note the input filename.
   std::cout << "\nimporting: " << ioparm.filename << std::endl;
*/

  // Builds filename strings.
  buildFilenames();

  // Gets adjacency matrix data.
  readInputFile<Graph>(ga, e, elist, e_wgt, unit);

  // ----------------------------------------------------------------------
  // Allocate data structures for main loop.
  num_nodes = new double [gparm.maxid];
  a = new double [gparm.maxid];
  Q = new double [gparm.n + 1];
  SV = new double [gparm.n + 1];
  joins = new apair  [gparm.n + 1];

  for (int i = 0; i < gparm.maxid; ++i)
  {
    a[i] = 0.0;
    num_nodes[i] = 1.0;
  }

  for (int i = 0; i < gparm.n + 1; ++i)
  {
    Q[i] = 0.0;
    joins[i].x = 0;
    joins[i].y = 0;
  }

  int t = 1;

  Qmax.y = -4294967296.0;
  Qmax.x = 0;
  SVmin.y = 4294967296.0;
  SVmin.x = 0;

#ifdef DEBUG
  std::cout << "now building initial dQ[]" << std::endl;
#endif

  // Builds dQ[] and h.
  buildDeltaQMatrix(cutstep, e, elist, e_wgt, unit);

  // Initialize f_joins, f_support files.
  if (WRITE_JOINS_FILE)
  {
    std::ofstream fjoins(ioparm.f_joins.c_str(), std::ios::trunc);
    fjoins << -1 << "\t" << -1 << "\t" << Q[0] << "\t0\n";
    fjoins.close();
  }

  if (ioparm.suppFlag)
  {
    std::ofstream fsupp(ioparm.f_support.c_str(), std::ios::trunc);
    dqSupport();
    fsupp << 0 << "\t" << supportTot << "\t" << supportAve << "\t" << 0
          << "\t->\t" << 0 << "\n";
    fsupp.close();
  }

  // ----------------------------------------------------------------------
  // Start FastCommunity algorithm.
#ifdef DEBUG
  std::cout << "starting algorithm now." << std::endl;
#endif

  size_type ord = num_vertices(ga);

  int* merge_l = new int[ord];
  int* merge_r = new int[ord];
  int npairs = 0;

  if (lf)
  {
#ifdef DEBUG
    printf("using leader file\n");
#endif

    gen_merges(ord, lf, merge_l, merge_r, npairs);

    for (int k = 0; k < npairs; ++k)
    {
#ifdef DEBUG
      printf("doing a booting merge\n");
      fflush(stdout);
#endif

      int i = merge_l[k];
      int j = merge_r[k];

      if (i == j) continue;

      // ---------------------------------
      // Find dQ.
      tuple_o dQnew;
      dQnew.m = dq[i].v->findItem(j)->heap_ptr->m;
      dQnew.sv_diff = dq[i].v->findItem(j)->heap_ptr->sv_diff;
#ifdef MB
      dQnew.MB_obj = 0;
#endif

#ifdef DEBUG
      std::cout << "Q[" << t - 1 << "] = " << Q[t - 1];
      std::cout << "SV[" << t - 1 << "] = " << SV[t - 1];
#endif

      // ---------------------------------
      // Merge the chosen communities.
      mergeCommunities(i, j);          // Merge community i into community j.
      joins[t].x = i;                  // Record merge of i(x) into j(y).
      joins[t].y = j;

      // Record Q(t).
#ifdef MB
      Q[t] = dQnew.MB_obj;
#else
      Q[t] = dQnew.m + Q[t - 1];
#endif

      SV[t] = dQnew.sv_diff + SV[t - 1];

      // ---------------------------------
      // Record join to file.
      if (WRITE_JOINS_FILE)
      {
        // Open file for writing the next join.
        std::ofstream fjoins(ioparm.f_joins.c_str(), std::ios::app);

        // Convert to external format.
        fjoins << joins[t].x - 1 << "\t" << joins[t].y - 1 << "\t";

        if ((Q[t] > 0.0 && Q[t] < 0.0000000000001) ||
            (Q[t] < 0.0 && Q[t] > -0.0000000000001))
        {
          fjoins << 0.0;
        }
        else
        {
          fjoins << Q[t];
        }

        fjoins << "\t" << t << "\n";
        fjoins.close();
      }

      // Note that it is the .joins file which contains both the dendrogram
      // and the corresponding Q values. The file format is tab-delimited
      // columns of data, where the columns are:
      //   1. The community which grows.
      //   2. The community which was absorbed.
      //   3. The modularity value Q after the join.
      //   4. The time step value.

#ifdef MB
      if (Q[t] < MB_THRESH)
      {
        Qmax.y = Q[t];
        Qmax.x = t;
      }
#else
      if (Q[t] > Qmax.y)
      {
        Qmax.y = Q[t];
        Qmax.x = t;
      }
#endif

      if (SV[t] < SVmin.y)
      {
        SVmin.y = SV[t];
        SVmin.x = t;
      }

      ++t;
    }
  }

  // ----------------------------------------------------------------------
  // Record some results.
  t2 = time(&t2);
  std::ofstream fout1(ioparm.f_parm.c_str(), std::ios::app);
  fout1 << "RESTART---:\t" << asctime(localtime(&t2));
  fout1.close();

#ifdef DEBUG
  printf("heapSize before start: %d\n", h->heapSize());
#endif

  int max_merges = h->heapSize();

  while (h->heapSize() >= 1)
  {
#ifdef DEBUG
    size_type num_merges = order - h->heapSize();

    if (num_merges % 10000 == 0)
    {
      printf("CNM merges: %lu\n", num_merges);
      fflush(stdout);
    }

    // ---------------------------------
    // Find largest dQ.
    if (ioparm.textFlag > 0)
    {
      h->printHeapTop10();

      std::cout << std::endl;

      for (int i = 0; i < gparm.maxid - 2; ++i)
      {
        if (dq[i].v)
        {
          printf("--------------------\n");
          dq[i].v->printHeap();
        }
      }
    }
#endif

    // Select maximum dQ_ij.  Convention: insert i into j.
    tuple_o dQmax = h->popMaximum();

    // No more joins possible.
    if (dQmax.m < -4000000000.0) break;

    if (dq[dQmax.i].v == NULL || dq[dQmax.j].v == NULL)
    {
      if (cutstep)
      {
        std::cout << "WARNING: invalid join (" << dQmax.i
                  << " " << dQmax.j << ") found at top of heap\n";
        continue;
      }
    }

    int isupport = dq[dQmax.i].v->returnNodecount();
    int jsupport = dq[dQmax.j].v->returnNodecount();

    // JWB: "Nodecount" is really neighboring adjacency count.  Hoping to
    //      make WT more effective, we looked instead at the number of nodes
    //      in the respective communities.
    //
    //      Tried: if (num_nodes[dQmax.i] < num_nodes[dQmax.j]).  Didn't help.
    if (isupport < jsupport)
    {
      // Merge community i into community j.
      mergeCommunities(dQmax.i, dQmax.j);

      // Record merge of i(x) into j(y).
      joins[t].x = dQmax.i;
      joins[t].y = dQmax.j;
    }
    else
    {
      // Take community j's heap pointer.
      dq[dQmax.i].heap_ptr = dq[dQmax.j].heap_ptr;

      // Mark it as i's.
      dq[dQmax.i].heap_ptr->i = dQmax.i;
      dq[dQmax.i].heap_ptr->j = dQmax.j;

      // Merge community j into community i.
      mergeCommunities(dQmax.j, dQmax.i);

      // Record merge of j(x) into i(y).
      joins[t].x = dQmax.j;
      joins[t].y = dQmax.i;
    }

    // Record Q(t).
#ifdef MB
    Q[t] = dQmax.MB_obj;
#else
    Q[t] = dQmax.m + Q[t - 1];
#endif

    // ---------------------------------
    // Record join to file.
    if (cutstep && WRITE_JOINS_FILE)
    {
      // Open file for writing the next join.
      std::ofstream fjoins(ioparm.f_joins.c_str(), std::ios::app);

      // Convert to external format.
      fjoins << joins[t].x - 1 << "\t" << joins[t].y - 1 << "\t";

      if ((Q[t] > 0.0 && Q[t] < 0.0000000000001) ||
          (Q[t] < 0.0 && Q[t] > -0.0000000000001))
      {
        fjoins << 0.0;
      }
      else
      {
        fjoins << Q[t];
      }

      fjoins << "\t" << t << " " << num_nodes[dQmax.i] << " "
             << num_nodes[dQmax.j] << std::endl;
      fjoins.close();
    }

    // Note that it is the .joins file which contains both the dendrogram and
    // the corresponding Q values. The file format is tab-delimited columns of
    // data, where the columns are:
    //   1. The community which grows.
    //   2. The community which was absorbed.
    //   3. The modularity value Q after the join.
    //   4. The time step value.

    if (t == ioparm.cutstep)
    {
      recordNetwork();
      recordGroupLists();
      groupListsStats();
    }

    // ---------------------------------
    // Record the support data to file.
    if (ioparm.suppFlag)
    {
      dqSupport();
      std::ofstream fsupp(ioparm.f_support.c_str(), std::ios::app);

      // time   remaining support   mean support   support_i --   support_j
      fsupp << t << "\t" << supportTot << "\t" << supportAve << "\t"
            << isupport;

      if (isupport < jsupport)
      {
        fsupp << "\t->\t";
      }
      else
      {
        fsupp << "\t<-\t";
      }

      fsupp << jsupport << "\n";
      fsupp.close();
    }

#ifdef MB
    if (Q[t] > MB_THRESH)
    {
      Qmax.y = Q[t];
      Qmax.x = t;
    }
#else
    if (Q[t] > Qmax.y)
    {
      Qmax.y = Q[t];
      Qmax.x = t;
    }
#endif

    ++t;
  }

  // ----------------------------------------------------------------------
  // Record some results.
  t1 = time(&t1);
  std::ofstream fout(ioparm.f_parm.c_str(), std::ios::app);
  fout << "---MODULARITY---\n";
  fout << "MAXQ------:\t" << Qmax.y << "\n";
  fout << "STEP------:\t" << Qmax.x << "\n";
  fout << "MINSV-----:\t" << SVmin.y << "\n";
  fout << "STEP------:\t" << SVmin.x << "\n";
  fout << "EXIT------:\t" << asctime(localtime(&t1));
  fout.close();

  int n = gparm.maxid - 2;
  if (gparm.n > n) n = gparm.n;

  for (int i = 0; i < n; ++i) leader[i] = i;

  // A value of 0 for num_communities indicates "let the algorithm decide the
  // stopping point."  Otherwise, stop at the number of communities given
  // by the user.
  int num_merges = num_communities == 0 ? Qmax.x + 1 :
                                          order - num_communities + 1;

  // Make sure the number of communities requested by the user doesn't cause
  // the algorithm to exceed the maximum number of merges.  The maximum number
  // of merges depends on the order and the number of components in the graph.
  if (num_merges > max_merges) num_merges = max_merges;
  num_communities = order - num_merges + 1;

  for (int i = 1; i < num_merges; ++i) leader[joins[i].x - 1] = joins[i].y - 1;

  delete [] merge_r;
  delete [] merge_l;
  delete [] e;
  delete [] elist;
  delete [] joins;
  delete [] a;
  delete [] num_nodes;
  delete [] Q;
  delete [] SV;

  for (int i = 0; i < gparm.maxid; ++i) if (dq[i].v) delete dq[i].v;
  delete [] dq;

  return Qmax.y;
}

template <typename etyp>
void buildDeltaQMatrix(bool cutstep, Edge<etyp>* e, Edge<etyp>* elist,
                       double* e_wgt, etyp unit)
{
  // Given that we've now populated a sparse (unordered) adjacency matrix e
  // (e), we now need to construct the intial dQ matrix according to the
  // definition of dQ which may be derived from the definition of modularity
  // Q:
  //
  //     Q(t) = \sum_{i} (e_{ii} - a_{i}^2) = Tr(e) - ||e^2||
  //
  // thus dQ is
  //
  //     dQ_{i,j} = 2* ( e_{i,j} - a_{i}a_{j} )
  //
  // where a_{i} = \sum_{j} e_{i,j} (i.e., the sum over the ith row).
  //
  // To create dQ, we must insert each value of dQ_{i,j} into a binary search
  // tree, for the jth column. That is, dQ is simply an array of such binary
  // search trees, each of which represents the dQ_{x,j} adjacency vector.
  // Having created dQ as such, we may happily delete the matrix e in order to
  // free up more memory.  The next step is to create a max-heap data
  // structure, which contains the entries of the following form
  // (value, s, t), where the heap-key is 'value'. Accessing the root of the
  // heap gives us the next dQ value, and the indices (s,t) of the vectors in
  // dQ which need to be updated as a result of the merge.


  // First we compute e_{i,j}, and the compute+store the a_{i} values. These
  // will be used shortly when we compute each dQ_{i,j}.
  Edge<etyp>* current;

  // Intially each e_{i,j} = 1/m.
//  double eij = (double) (0.5 / gparm.m);

  // JWB: In the weighted case, e_{i,j} = w_ij/2T, so the declaration above
  //      is no longer valid.  We'll amend below.

#ifdef DEBUG
  printf("gp.T: %f\n", gparm.T);
#endif

  // For each row.
  for (int i = 1; i < gparm.maxid; ++i)
  {
    a[i] = 0.0;

    // Ensure it exists.
    if (e[i].so != 0)
    {
      // Grab first edge.
      current = &e[i];

      // JWB: Now, find the weighted e_ij (recall 0-based).
      double w_ij = e_wgt[current->eid];
      a[i] += w_ij / (2 * gparm.T);

      // Loop through remaining edges.
      while (current->next != NULL)
      {
        current = current->next;
        double w_ij = e_wgt[current->eid];
        a[i] += w_ij / (2 * gparm.T);
      }

      // Calculate initial value of Q.
      Q[0] += -1.0 * a[i] * a[i];
    }
  }

  // Now, we create an empty (ordered) sparse matrix dq[].
  dq = new nodenub [gparm.maxid];

  for (int i = 0; i < gparm.maxid; ++i)
  {
    // No pointer in the heap at first.
    dq[i].heap_ptr = NULL;

    if (e[i].so != 0)
    {
      dq[i].v = new vektor(i, 2 + (int) floor(gparm.m * a[i]));
    }
    else
    {
      dq[i].v = NULL;
    }
  }

  // Allocate max-heap of size = number of nodes.
  h = new maxheap(gparm.n);

  // Now, we do all the work, which happens as we compute and insert each
  // dQ_{i,j} into the corresponding (ordered) sparse vector dq[i]. While
  // computing each dQ for a row i, we track the maximum dQmax of the row and
  // its (row,col) indices (i,j). Upon finishing all dQ's for a row, we
  // insert the tuple into the max-heap hQmax. That insertion returns the
  // itemaddress, which we then store in the nodenub heap_ptr for that row's
  // vector.
  double dQ;
  tuple_o dQmax;         // For heaping the row maxes.

  for (int i = 1; i < gparm.maxid; ++i)
  {
    if (e[i].so != 0)
    {
      // Grab first edge.
      current = &e[i];
      double w_ij = e_wgt[current->eid];

      dQ = 2 * (w_ij / (2 * gparm.T) - (a[current->so] * a[current->si]));

      dQmax.m = dQ;                          // Assume it is maximum so far.
      dQmax.i = current->so;                 // Store its (row,col).
      dQmax.j = current->si;
      dq[i].v->insertItem(current->si, dQ);  // Insert its dQ.

      while (current->next != NULL)
      {
        // Step to next edge.
        current = current->next;
        double w_ij = e_wgt[current->eid];

        dQ = 2 * (w_ij / (2 * gparm.T) - (a[current->so] * a[current->si]));

        // If dQ larger than current max.
        if (dQ > dQmax.m)
        {
          dQmax.m = dQ;            // Replace it as maximum so far,
          dQmax.j = current->si;   // and store its (col).
        }

        // Insert it into vector[i]
        dq[i].v->insertItem(current->si, dQ);
      }

      // Store the pointer to its loc in heap.
      dq[i].heap_ptr = h->insertItem(dQmax);
    }
  }

  // Free-up adjacency matrix memory in two shots.
  if (!cutstep) delete [] elist;
  if (!cutstep) delete [] e;

  return;
}

void buildFilenames()
{
  ioparm.f_input = ioparm.d_in + ioparm.filename;
  ioparm.f_parm = ioparm.d_out + ioparm.s_scratch + "-fc_" +
                  ioparm.s_label + ".info";
  ioparm.f_joins = ioparm.d_out + ioparm.s_scratch + "-fc_" +
                   ioparm.s_label + ".joins";
  ioparm.f_support = ioparm.d_out + ioparm.s_scratch + "-fc_" +
                     ioparm.s_label + ".supp";
  ioparm.f_net = ioparm.d_out + ioparm.s_scratch + "-fc_" +
                 ioparm.s_label + ".wpairs";
  ioparm.f_group = ioparm.d_out + ioparm.s_scratch + "-fc_" +
                   ioparm.s_label + ".groups";
  ioparm.f_gstats = ioparm.d_out + ioparm.s_scratch + "-fc_" +
                    ioparm.s_label + ".hist";

  if (true)
  {
    std::ofstream flog(ioparm.f_parm.c_str(), std::ios::trunc);
    flog.close();
  }

  time_t t;
  t = time(&t);

  std::ofstream flog(ioparm.f_parm.c_str(), std::ios::app);

  flog << "FASTCOMMUNITY_INFERENCE_ALGORITHM\n";
  flog << "START-----:\t" << asctime(localtime(&t));
  flog << "---FILES--------\n";
  flog << "DIRECTORY-:\t" << ioparm.d_out << "\n";
  flog << "F_IN------:\t" << ioparm.filename << "\n";
  flog << "F_JOINS---:\t"
       << ioparm.s_scratch + "-fc_" + ioparm.s_label + ".joins" << "\n";
  flog << "F_INFO----:\t"
       << ioparm.s_scratch + "-fc_" + ioparm.s_label + ".info" << "\n";

  if (ioparm.suppFlag)
  {
    flog << "F_SUPP----:\t"
         << ioparm.s_scratch + "-fc_" + ioparm.s_label + ".supp" << "\n";
  }

  if (ioparm.cutstep > 0)
  {
    flog << "F_NET-----:\t"
         << ioparm.s_scratch + "-fc_" + ioparm.s_label + ".wpairs" << "\n";
    flog << "F_GROUPS--:\t"
         << ioparm.s_scratch + "-fc_" + ioparm.s_label + ".groups" << "\n";
    flog << "F_GDIST---:\t"
         << ioparm.s_scratch + "-fc_" + ioparm.s_label + ".hist" << "\n";
  }

  flog.close();

  return;
}

// ---------------------------------------------------------------------------
// Returns the support of the dQ[].
void dqSupport()
{
  int total = 0;
  int count = 0;

  for (int i = 0; i < gparm.maxid; ++i)
  {
    if (dq[i].heap_ptr != NULL)
    {
      total += dq[i].v->returnNodecount();
      ++count;
    }
  }

  supportTot = total;
  supportAve = total / (double) count;

  return;
}

template <typename etyp>
void groupListsSetup(Edge<etyp>* e)
{
  c = new stub [gparm.maxid];

  for (int i = 0; i < gparm.maxid; ++i)
  {
    // Note: internal indexing.
    if (e[i].so != 0)
    {
      cnmlist* newList = new cnmlist;   // Create new community member,
      newList->index = i;               //   with index i.

      c[i].members = newList;           // Point ith community at newList.
      c[i].size = 1;                    // Point ith community at newList.
      c[i].last = newList;              // Point last[] at that element, too.
      c[i].valid = true;                // Mark as valid community.
      c[i].one = true;
    }
  }

  return;
}

// ---------------------------------------------------------------------------
// Function for computing statistics on the list of groups.
void groupListsStats()
{
  gstats.numgroups = 0;
  gstats.maxsize = 0;
  gstats.minsize = gparm.maxid;

  double count = 0.0;

  for (int i = 0; i < gparm.maxid; ++i)
  {
    if (c[i].valid)
    {
      // Count number of communities.
      ++gstats.numgroups;
      count += 1.0;

      // Find biggest community.
      if (c[i].size > gstats.maxsize) gstats.maxsize = c[i].size;

      // Find smallest community.
      if (c[i].size < gstats.minsize) gstats.minsize = c[i].size;

      // Compute mean group size.
      gstats.meansize = (double) (c[i].size) / count +
                        (((double) (count - 1.0) / count) * gstats.meansize);
    }
  }

  gstats.sizehist = new double [gstats.maxsize + 1];

  for (int i = 0; i < gstats.maxsize + 1; ++i) gstats.sizehist[i] = 0;

  count = 0.0;

  for (int i = 0; i < gparm.maxid; ++i)
  {
    if (c[i].valid)
    {
      // Tabulate histogram of sizes.
      gstats.sizehist[c[i].size] += 1.0;
      count += 1.0;
    }
  }

  // Convert histogram to pdf, and write it to disk.
  for (int i = 0; i < gstats.maxsize + 1; ++i)
  {
    gstats.sizehist[i] = gstats.sizehist[i] / count;
  }

  std::ofstream fgstat(ioparm.f_gstats.c_str(), std::ios::trunc);
  for (int i = gstats.minsize; i < gstats.maxsize + 1; ++i)
  {
    fgstat << i << "\t" << gstats.sizehist[i] << "\n";
  }
  fgstat.close();

  // Record some statistics.
  time_t t1;
  t1 = time(&t1);

  std::ofstream fout(ioparm.f_parm.c_str(), std::ios::app);
  fout << "---GROUPS-------\n";
  fout << "NUMGROUPS-:\t" << gstats.numgroups << "\n";
  fout << "MINSIZE---:\t" << gstats.minsize << "\n";
  fout << "MEANSIZE--:\t" << gstats.meansize << "\n";
  fout << "MAXSIZE---:\t" << gstats.maxsize << "\n";
  fout.close();

  return;
}

void groupListsUpdate(const int x, const int y)
{
  c[y].last->next = c[x].members;   // Attach c[y] to end of c[x].
  c[y].last = c[x].last;            // Update last[] for community y.
  c[y].size += c[x].size;           // Add size of x to size of y.

  c[x].members = NULL;              // Delete community[x].
  c[x].valid = false;
  c[x].size = 0;
  c[x].last = NULL;                 // Delete last[] for community x.

  return;
}

void mergeCommunities(int i, int j)
{
  std::list<int> xvals;

  num_nodes[j] += num_nodes[i];

  // Go through i's neighbors. If we find one that isn't adjacent
  // to j, we know the dQ by calculation.
  dpair* list = dq[i].v->returnTreeAsList();
  dpair* current = list;

  while (current != NULL)
  {
    int x = current->x;
    if (x != j)
    {
      xvals.push_back(x);
      if (!dq[j].v->findItem(x))
      {
        double dqjx = -2 * a[j] * a[x];
        dq[j].v->insertItem(x, dqjx);
        dq[x].v->insertItem(j, dqjx);
      }
      dq[x].v->deleteItem(i);
    }

    dpair* temp = current;
    current = current->next;
    delete temp;
  }

  // Go through j's neighbors. If we find one that isn't adjacent
  // to i, we know the dQ by calculation.
  list = dq[j].v->returnTreeAsList();
  current = list;

  while (current != NULL)
  {
    int x = current->x;

    if (x != i)
    {
      xvals.push_back(x);

      if (!dq[i].v->findItem(x))
      {
        double dqix = -2 * a[i] * a[x];
        dq[i].v->insertItem(x, dqix);
      }
    }

    dpair* temp = current;
    current = current->next;
    delete temp;
  }

  a[j] += a[i];

  xvals.sort();
  xvals.unique();

  // Go through all neighbors.
  // dQ_x'[j] = dQ'_j[x] = dQ_i[x] + dQ_j[x]
  // There is no need to break out JIX, IJX, etc.
  std::list<int>::iterator lit = xvals.begin();
  for ( ; lit != xvals.end(); ++lit)
  {
    int x = *lit;
    double dqix = dq[i].v->findItem(x)->heap_ptr->m;
    double dqjx = dq[j].v->findItem(x)->heap_ptr->m;

    dq[j].v->insertItem(x, dqix + dqjx);
    dq[x].v->insertItem(j, dqix + dqjx);

    tuple_o newMax = dq[x].v->returnMaxStored();
    h->updateItem(dq[x].heap_ptr, newMax);
  }

  dq[j].v->deleteItem(i);
  tuple_o newMax = dq[j].v->returnMaxStored();
  h->updateItem(dq[j].heap_ptr, newMax);

  delete dq[i].v;

  dq[i].v = NULL;
  dq[i].heap_ptr = NULL;
}

bool parseCommandLine(int argc, char* argv[])
{
  // Description of commandline arguments
  // -f <filename>   Give the target .pairs file to be processed
  // -l <text>       The text label for this run; used to build output
  //                 filenames.
  // -t <int>        Period of timer for reporting progress of computation
  //                 to screen.
  // -s              Calculate and track the support of the dQ matrix.
  // -v --v ---v     Differing levels of screen output verbosity.
  // -c <int>        Record the aglomerated network at step <int>.

  // If no arguments, return statement about program usage.
  if (argc <= 1)
  {
    std::cout << "\nThis program runs the fast community structure inference "
                 "algorithm due to Clauset, Newman and Moore on an input "
                 "graph in the .pairs format. This version is the full "
                 "max-heap version originally described in cond-mat/0408187. "
                 "The program requires the input network connectivity to be "
                 "formatted in the following specific way: the graph must be "
                 "simple and connected, where each edge is written on a line "
                 "in the format 'u v' (e.g., 3481  3483).\n";

    std::cout << "To run the program, you must specify the input network file "
                 "(-f file.pairs). Additionally, you can differentiate runs "
                 "on the same input file with a label l test_run) which is "
                 "imbedded in all corresponding output files. Because the "
                 "algorithm is deterministic, you can specify a point (-c C) "
                 "at which to cut the dendrogram; the program will write out "
                 "various information about the clustered network: a list of "
                 "clusters, the clustered connectivity, and the cluster size "
                 "distribution. Typically, one wants this value to be the "
                 " time at which modularity Q was maximized (that time is "
                 "recorded in the .info file for a given run).\n";

    std::cout << "Examples:\n"
                 "  ./FastCommunity -f network.pairs -l test_run\n"
                 "  ./FastCommunity -f network.pairs -l test_run -c 1232997\n"
                 "\n";

    return false;
  }

  int argct = 1;

  while (argct < argc)
  {
    std::string temp = argv[argct];

    if (temp == "-files")
    {
      std::cout << "\nBasic files generated:\n"
                   "-- .INFO\n"
                   "   Various information about the program's running. "
                   "Includes a listing of the files it generates, number of "
                   "vertices and edges processed, the maximum modularity "
                   "found and the corresponding step (you can re-run the "
                   "program with this value in the -c argument to have it "
                   "output the contents of the clusters, etc. when it "
                   "reaches that step again (not the most efficient "
                   "solution, but it works)), start/stop time, and when -c "
                   "is used, it records some information about the "
                   "distribution of cluster sizes.\n";

      std::cout << "-- .JOINS\n"
                   "   The dendrogram and modularity information from the "
                   "algorithm. The file format is tab-delimited columns of "
                   "data, where the columns are:\n"
                   " 1. the community index which absorbs\n"
                   " 2. the community index which was absorbed\n"
                   " 3. the modularity value Q after the join\n"
                   " 4. the time step of the join\n";

      std::cout << "\nOptional files generated (at time t=C when -c C "
                   "argument used):\n"
                   "-- .WPAIRS\n"
                   "   The connectivity of the clustered graph in a .wpairs "
                   "file format (i.e., weighted edges). The edge weights "
                   "should be the dQ values associated with that clustered "
                   "edge at time C. From this format, it's easy to convert "
                   "into another for visualization (e.g., pajek's .net "
                   "format).\n";

      std::cout << "-- .HIST\n"
                   "   The size distribution of the clusters.\n";

      std::cout << "-- .GROUPS\n"
                   "   A list of each group and the names of the vertices "
                   "which compose it (this is particularly useful for "
                   "verifying that the clustering makes sense - tedious but "
                   "important).\n\n";

      return false;
    }
    else if (temp == "--filename")            // Input file name.
    {
      ++argct;
      temp = argv[argct];
      std::string ext = ".g";

      std::string::size_type pos = temp.find(ext, 0);
      if (pos == std::string::npos)
      {
        std::cout << " Error: Input file must have terminating .g extension.\n";
        return false;
      }

      ext = "/";
      pos = std::string::npos;
      for (std::string::size_type i = 0; i < temp.size(); ++i)
      {
        if (temp[i] == '/') pos = i;
      }

      if (pos == std::string::npos)
      {
        ioparm.d_in = "";
        ioparm.filename = temp;
      }
      else
      {
        ioparm.d_in = temp.substr(0, pos + 1);
        ioparm.filename = temp.substr(pos + 1, temp.size() - pos - 1);
      }

      ioparm.d_out = ioparm.d_in;

      // Now grab the filename sans extension for building outputs files.
      for (std::string::size_type i = 0; i < ioparm.filename.size(); ++i)
      {
        if (ioparm.filename[i] == '.') pos = i;
      }

      ioparm.s_scratch = ioparm.filename.substr(0, pos);
    }
    else if (temp == "-l")            // s_label.
    {
      ++argct;

      if (argct < argc)
      {
        ioparm.s_label = argv[argct];
      }
      else
      {
        std::cout << " Warning: missing modifier for -l argument; using "
                     "default.\n";
      }
    }
    else if (temp == "-t")            // Timer value.
    {
      ++argct;

      if (argct < argc)
      {
        ioparm.timer = atoi(argv[argct]);

        std::cout << ioparm.timer << std::endl;

        if (ioparm.timer == 0 || strlen(argv[argct]) > temp.length())
        {
          std::cout << " Warning: malformed modifier for -t; using default.\n";
          --argct;
          ioparm.timer = 20;
        }
      }
      else
      {
        std::cout << " Warning: missing modifier for -t argument; using "
                     "default.\n";
        --argct;
      }
    }
    else if (temp == "-c")            // cut value
    {
      ++argct;

      if (argct < argc)
      {
        ioparm.cutstep = atoi(argv[argct]);

        if (ioparm.cutstep == 0)
        {
          std::cout << " Warning: malformed modifier for -c; disabling "
                       "output.\n";
          --argct;
        }

        ++ioparm.cutstep;
      }
      else
      {
        std::cout << " Warning: missing modifier for -t argument; using "
                     "default.\n";
        --argct;
      }
    }
    else if (temp == "-s")
    {
      ioparm.suppFlag = true;
    }
    else if (temp == "-v")
    {
      ioparm.textFlag = 1;
    }
    else if (temp == "--v")
    {
      ioparm.textFlag = 2;
    }
    else if (temp == "---v")
    {
      ioparm.textFlag = 3;
    }
    else
    {
      std::cout << "Unknown commandline argument: " << argv[argct]
                << std::endl;
    }

    ++argct;
  }

  return true;
}

// ---------------------------------------------------------------------------
// Records the agglomerated list of indices for each valid community.
void recordGroupLists()
{
  std::ofstream fgroup(ioparm.f_group.c_str(), std::ios::trunc);

  for (int i = 0; i < gparm.maxid; ++i)
  {
    if (c[i].valid)
    {
      // External format.
      fgroup << "GROUP[ " << i - 1 << " ][ " << c[i].size << " ]\n";

      cnmlist* current = c[i].members;
      while (current != NULL)
      {
        // External format.
        fgroup << current->index - 1 << "\n";

        current = current->next;
      }
    }
  }

  fgroup.close();

  return;
}

// ---------------------------------------------------------------------------
// Records the network as currently agglomerated.
void recordNetwork()
{
  std::ofstream fnet(ioparm.f_net.c_str(), std::ios::trunc);

  for (int i = 0; i < gparm.maxid; ++i)
  {
    if (dq[i].heap_ptr != NULL)
    {
      // Get a list of items in dq[i].v.
      dpair* list = dq[i].v->returnTreeAsList();

      // Store ptr to head of list.
      dpair* current = list;
      while (current != NULL)
      {
        // source        target        weight    (external representation)
        fnet << i - 1 << "\t" << current->x - 1 << "\t" << current->y << "\n";

        // Clean up memory and move to next.
        dpair* temp = current;
        current = current->next;
        delete temp;
      }
    }
  }

  fnet.close();

  return;
}

// ***********************************************************************
// gen_merges()
//
// Given a file of "id leader" relationships, compute a sequence of
// merges that establish those relationships.  Order matters; a vertex
// can appear on the left-hand side of a merge only once -- then it is
// merged into its leader and no longer has an identity in cnm.
//
// This is used by the -lf option, which reads an "id leader" file.
// ***********************************************************************
void gen_merges(int nodes, int* id, int* leader, int& npairs)
{
  int* result_id = new int[nodes];
  int* result_leader = new int[nodes];

  printf("gen_merges(%d,..,..,%d)\n", nodes, npairs);
  fflush(stdout);

  int* aleader = (int*) malloc(nodes * sizeof(int));
  int* parent = (int*) malloc(nodes * sizeof(int));

  for (int i = 0; i < nodes; ++i)
  {
    aleader[i] = 0;
    parent[i] = 0;
  }

  int count = 0;

  for (int i = 0; i < npairs; ++i)
  {
    ++aleader[leader[i] - 1];
    ++count;
    parent[id[i] - 1] = leader[i];
  }

  int ind = 0;

  while (count)
  {
    for (int i = 0; i < nodes; ++i)
    {
      if (parent[i] && !aleader[i])
      {
        printf("%d %d\n", i, parent[i]);
        result_id[ind] = i;
        result_leader[ind] = parent[i];
        --aleader[parent[i] - 1];
        --count;
        parent[i] = 0;
        ++ind;
      }
    }
  }

  for (int i = 0; i < ind; ++i)
  {
    id[i] = result_id[i] + 1;
    leader[i] = result_leader[i];
    printf("merge: (%d,%d)\n", id[i], leader[i]);
    fflush(stdout);
  }

  npairs = ind;

  free(aleader);
  free(parent);
  delete [] result_id;
  delete [] result_leader;
}

void gen_merges(int nodes, char* fname, int* id, int* leader, int& npairs)
{
  FILE* file;

  if ((file = fopen(fname, "r")) == NULL)
  {
    printf("Cannot open file %s", fname);
    return;
  }

  npairs = 0;
  while (2 == fscanf(file, "%d%d", &id[npairs], &leader[npairs])) ++npairs;

  gen_merges(nodes, id, leader, npairs);
}

void read_merges(int nodes, char* fname, int* id, int* leader, int& npairs)
{
  FILE* file;

  if ((file = fopen(fname, "r")) == NULL)
  {
    printf("Cannot open file %s", fname);
    return;
  }

  npairs = 0;
  while (2 == fscanf(file, "%d%d", &id[npairs], &leader[npairs])) ++npairs;
}

template <typename Graph>
double
wcnm(Graph& ga, typename graph_traits<Graph>::size_type* leader,
     int& num_communities, int iter = 1, double* e_weight = 0)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;

  typedef xmt_hash_table<int64_t, size_type> hash_table_t;
  typedef xmt_hash_table<int64_t, double> dhash_table_t;

  size_type order = num_vertices(ga);
  size_type size = num_edges(ga);

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, ga);
  edge_id_map<Graph> eid_map = get(_edge_id_map, ga);

#ifdef DEBUG
  printf("order: %lu\n", order);
  printf(" size: %lu\n", size);
#endif

  // There is a maximum of order communities.  Make sure num_communities
  // isn't larger than this.
  if (num_communities > static_cast<int>(order)) num_communities = order;

  bool i_own_e_weight = (e_weight == 0);
  if (i_own_e_weight) e_weight = (double*) malloc(sizeof(double) * size);

//  dhash_table_t e_weight((int)(ga.get_size() * 1.5));

  edge_iterator eiter = edges(ga);

  if (i_own_e_weight)
  {
    for (size_type i = 0; i < size; ++i)
    {
      size_type eid = get(eid_map, eiter[i]);
      e_weight[eid] = 1.0;
    }
  }

  mt_timer timer;
  weight_by_neighborhoods<Graph> wbn(ga, e_weight);

  timer.start();
  wbn.run(iter);
  timer.stop();

#ifdef DEBUG
  printf("edge weighting time: %f\n", timer.getElapsedSeconds());
#endif

#if WRITE_WEIGHTED_DIMACS
  std::ofstream fwdimacs("weighted.dimacs", std::ios::trunc);
  fwdimacs << "p sp " << order << " " << size << std::endl;
#endif

  W = 0;

  for (size_type i = 0; i < size; ++i)
  {
    edge_descriptor e = eiter[i];
    size_type eid = get(eid_map, e);
    double wgt = e_weight[eid];

    if (wgt < 0.0)
    {
      printf("edge(%lu): bad weight: %f\n", eid, wgt);
      exit(1);
    }

    W += wgt;

#if WRITE_WEIGHTED_DIMACS
    vertex_descriptor v1 = source(e, ga);
    vertex_descriptor v2 = target(e, ga);
    size_type v1id = get(vid_map, v1);
    size_type v2id = get(vid_map, v2);

    printf("wcnm wgt (%lu,%lu): %f\n", v1id, v2id, wgt);
    double scaled_wgt = wgt * 1e7;
    fwdimacs << "a " << v1id + 1 << " " << v2id + 1 << " " << (int) scaled_wgt
             << std::endl;
#endif
  }

#ifdef DEBUG
  printf("W: %f\n", W);
#endif

  size_type* parent = new size_type[order];
  for (size_type i = 0; i < order; ++i) parent[i] = i;

  double modularity = fastcommunity_mh(ga, e_weight, parent, 0,
                                       num_communities, 0, true);

  for (size_type i = 0; i < order; ++i)
  {
    size_type p = i;
    while (p != parent[p]) p = parent[p];

    leader[i] = p;
  }

  delete [] parent;

  if (i_own_e_weight) free(e_weight);

  return modularity;
}

}

#undef MB_THRESH
#undef WRITE_JOINS_FILE
#undef WRITE_WEIGHTED_DIMACS

#endif
