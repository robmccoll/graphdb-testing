#include  <cstdio>
#include  <omp.h>
#include  <vector>
#include <sys/time.h>
#include <sys/resource.h>

#include    "mtgl/adjacency_list.hpp"
#include    "mtgl/breadth_first_search.hpp"
#include    "mtgl/pagerank.hpp"

extern "C" {
#include  "timer.h"
}

#define E_A(X,...) fprintf(stderr, "%s %s %d:\n\t" #X "\n", __FILE__, __func__, __LINE__, __VA_ARGS__); exit(-1);
#define E(X) E_A(X,NULL)
#define V_A(X,...) fprintf(stdout, "%s %s %d:\n\t" #X "\n", __FILE__, __func__, __LINE__, __VA_ARGS__);
#define V(X) V_A(X,NULL)
#define R_A(X,...) fprintf(stdout, "RSLT: " X, __VA_ARGS__);
#define R(X) R_A(X,NULL)

using namespace mtgl;

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * REIMPLEMENTATION OF SHILOACH-VISHKIN 
 * I could not get their version to work with their adjacency_list 
 * (including compiling their test with the adj. switched onlist in the svn 
 * as of 2013/02/26 or download v1.1.*)
 * So I copied theirs, OpenMP'ed it, copied in their initialization, and added
 * the reduction of the number of components.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
template <typename Graph, typename ComponentMap>
typename graph_traits<Graph>::size_type
shiloach_vishkin(Graph& g, ComponentMap& result)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;
  typedef typename graph_traits<Graph>::thread_edge_iterator
          thread_edge_iterator;

  vertex_id_map<Graph> vid_map = get(_vertex_id_map, g);
  vertex_property_map<Graph, vertex_descriptor> D(g);

  size_type size = num_edges(g);
  size_type order = num_vertices(g);

  size_type num_streams;
  #pragma omp parallel 
  {
    #pragma omp master
    num_streams = omp_get_num_threads();
  }

  // Initialize each vertex to have itself as its leader.
  #pragma omp parallel
  {
    size_type stream_id = omp_get_thread_num();
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator verts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++verts)
    {
      put(result, *verts, start_pos);
    }
  }

  #pragma omp parallel
  {
    size_type stream_id = omp_get_thread_num();
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator verts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++verts)
    {
      vertex_descriptor v = *verts;
      D[v] = *(thread_vertices(result[v], g));
    }
  }

  int graft = 1;

  while (graft != 0)
  {
    graft = 0;

    #pragma omp parallel
    {
      size_type stream_id = omp_get_thread_num();
      size_type start_pos = begin_block_range(size, stream_id, num_streams);
      size_type end_pos = end_block_range(size, stream_id, num_streams);

      thread_edge_iterator edgs = thread_edges(start_pos, g);
      for ( ; start_pos != end_pos; ++start_pos, ++edgs)
      {
        edge_descriptor e = *edgs;

	vertex_descriptor u = source(e, g);
	vertex_descriptor v = target(e, g);

	if (get(vid_map, D[u]) < get(vid_map, D[v]) && D[v] == D[D[v]])
	{
	  D[D[v]] = D[u];
	  graft = 1;
	}

	if (get(vid_map, D[v]) < get(vid_map, D[u]) && D[u] == D[D[u]])
	{
	  D[D[u]] = D[v];
	  graft = 1;
	}
      }
    }

    #pragma omp parallel
    {
      size_type stream_id = omp_get_thread_num();
      size_type start_pos = begin_block_range(order, stream_id, num_streams);
      size_type end_pos = end_block_range(order, stream_id, num_streams);

      thread_vertex_iterator verts = thread_vertices(start_pos, g);
      for ( ; start_pos != end_pos; ++start_pos, ++verts)
      {
        vertex_descriptor v = *verts;

        while (D[v] != D[D[v]])
        {
          D[v] = D[D[v]];
        }
      }
    }
  }

  #pragma omp parallel
  {
    size_type stream_id = omp_get_thread_num();
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator verts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++verts)
    {
      vertex_descriptor v = *verts;
      result[v] = get(vid_map, D[v]);
    }
  }

  size_type count = 0;
  #pragma omp parallel reduction(+:count)
  {
    size_type stream_id = omp_get_thread_num();
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator verts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++verts)
    {
      count += (start_pos == get(result, *verts));
    }
  }

  return count;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * DISTANCE BFS VISITOR
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
template<typename Graph, typename DistMap>
class hop_count_bfs_visitor : public default_bfs_visitor<Graph> {
  public:
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

    hop_count_bfs_visitor(DistMap &d) : dist(d) { }

    void tree_edge(edge_descriptor & e, Graph & g) {
      vertex_descriptor src = source(e, g);
      vertex_descriptor dst = target(e, g);

      mt_write(dist[dst], dist[src] + 1);
    }

  private:
    DistMap & dist;
};

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * PAGERANK
 * Copy / reimplement with OpenMP
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
template <typename Graph, typename RankMap>
int pagerank_omp(Graph& g, RankMap& rank,
              double delta = .00001, double dampen = .8, int maxiter = 100)
{
  typedef typename graph_traits<Graph>::size_type size_type;
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::thread_vertex_iterator
          thread_vertex_iterator;

  size_type order = num_vertices(g);

  typedef vertex_property_map<Graph, double> AccMap;
  AccMap acc(g);

  size_type num_streams = 1;
  #pragma omp parallel 
  {
    #pragma omp master
    num_streams = omp_get_num_threads();
  }

  #pragma omp parallel
  {
    size_type stream_id = omp_get_thread_num();
    size_type start_pos = begin_block_range(order, stream_id, num_streams);
    size_type end_pos = end_block_range(order, stream_id, num_streams);

    thread_vertex_iterator tverts = thread_vertices(start_pos, g);
    for ( ; start_pos != end_pos; ++start_pos, ++tverts)
    {
      rank[*tverts] = 1.0 / ((double)order);
    }
  }

  // Get the accumulation array of the vertex degrees.
  //   undirectedS, directedS - out_degree()
  //   bidirectionalS - in_degree()
  size_type* accum_deg = (size_type*) malloc(sizeof(size_type) * (order + 1));
  detail::accumulate_dir_degree<Graph> add;
  add(accum_deg, g);

  size_type num_blocks = 1;

  int iter_cnt = 0;
  double maxdiff = 0.0;

  do
  {
    ++iter_cnt;

    double sum = 0.0;

    #pragma omp parallel reduction(+:sum)
    {
      size_type stream_id = omp_get_thread_num();
      size_type start_pos = begin_block_range(order, stream_id, num_streams);
      size_type end_pos = end_block_range(order, stream_id, num_streams);

      thread_vertex_iterator tverts = thread_vertices(start_pos, g);
      for ( ; start_pos != end_pos; ++start_pos, ++tverts)
      {
        vertex_descriptor v = *tverts;

        acc[v] = 0.0;
        sum += rank[v];
      }
    }

    detail::compute_acc<Graph, RankMap, AccMap> cacc;
    cacc(g, rank, acc, accum_deg, num_blocks);

    double adjustment = (1 - dampen) / order * sum;

    // Adjustment for zero-outdegree vertices.
    sum = 0;

    #pragma omp parallel reduction(+:sum)
    {
      size_type stream_id = omp_get_thread_num();
      size_type start_pos = begin_block_range(order, stream_id, num_streams);
      size_type end_pos = end_block_range(order, stream_id, num_streams);

      thread_vertex_iterator tverts = thread_vertices(start_pos, g);
      for ( ; start_pos != end_pos; ++start_pos, ++tverts)
      {
        vertex_descriptor v = *tverts;
        sum += (out_degree(v, g) == 0) * rank[v];
      }
    }

    adjustment += dampen * sum / order;

    // Compute new solution vector and scaling factor.
    double norm = 0.0;

    #pragma omp parallel
    {
      size_type stream_id = omp_get_thread_num();
      size_type start_pos = begin_block_range(order, stream_id, num_streams);
      size_type end_pos = end_block_range(order, stream_id, num_streams);

      thread_vertex_iterator tverts = thread_vertices(start_pos, g);
      for ( ; start_pos != end_pos; ++start_pos, ++tverts)
      {
        vertex_descriptor v = *tverts;

        acc[v] = adjustment + (dampen * acc[v]);

        double tmp = acc[v] >= 0 ? acc[v] : -1 * acc[v];
        if (tmp > norm) {
	  #pragma omp critical
	  if (tmp > norm) {
	    norm = tmp;
	  }
	}
      }
    }

    maxdiff = 0;

    #pragma omp parallel
    {
      size_type stream_id = omp_get_thread_num();
      size_type start_pos = begin_block_range(order, stream_id, num_streams);
      size_type end_pos = end_block_range(order, stream_id, num_streams);

      thread_vertex_iterator tverts = thread_vertices(start_pos, g);
      for ( ; start_pos != end_pos; ++start_pos, ++tverts)
      {
        vertex_descriptor v = *tverts;

        double oldval = rank[v];
        double newval = acc[v] / norm;
        rank[v] = newval;

        double absdiff = oldval > newval ? oldval - newval : newval - oldval;
        if (absdiff > maxdiff) {
	  #pragma omp critical
	  if(absdiff > maxdiff) {
	    maxdiff = absdiff;
	  }
	}
      }
    }

  } while (maxdiff > delta && iter_cnt < maxiter);

  free(accum_deg);
  return iter_cnt;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
 * MAIN
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

int main(int argc, char *argv[]) {
  if(argc < 3) {
    E_A(Not enough arguments. Usage %s graphfile actionsfile, argv[0]);
  }

  R("{\n")
  R("\"type\":\"mtgl\",\n")

  FILE * fp = fopen(argv[1], "r");

  int64_t nv;
  int64_t ne;
  int64_t * off, * ind, * wgt;

  const uint64_t endian_check = 0x1234ABCDul;
  uint64_t check;

  fread(&check, sizeof(uint64_t), 1, fp);

  if(check != endian_check) {
    E(Endianness does not agree.  Order swapping not implemented.);
  }

  fread(&nv, sizeof(int64_t), 1, fp);
  fread(&ne, sizeof(int64_t), 1, fp);

  off = (int64_t *)malloc(sizeof(int64_t) * nv+1);
  ind = (int64_t *)malloc(sizeof(int64_t) * ne);
  wgt = (int64_t *)malloc(sizeof(int64_t) * ne);

  fread(off, sizeof(int64_t), nv+1, fp);
  fread(ind, sizeof(int64_t), ne, fp);
  fread(wgt, sizeof(int64_t), ne, fp);

  fclose(fp);

  R_A("\"nv\":%ld,\n", nv)
  R_A("\"ne\":%ld,\n", ne)
  R("\"results\": {\n")

  V(Creating graph...);

  typedef adjacency_list<undirectedS> Graph;
  typedef graph_traits<Graph>::size_type size_type;
  typedef graph_traits<Graph>::vertex_iterator vertex_iterator;

  size_type * src = new size_type[ne];
  size_type * dst = new size_type[ne];

  for(uint64_t e = 0; e < ne; e++) {
    src[e] = 0;
    dst[e] = 0;
  }


  V(Loading data into graph...);
  tic();
  for(uint64_t v = 0; v < nv; v++) {
    for(uint64_t i = off[v]; i < off[v+1]; i++) {
      src[i] = v; dst[i] = ind[i];
    }
  }

  Graph g;
  init(nv, ne, src, dst, g);

  double build_time = toc();
  R("\"build\": {\n")
  R("\"name\":\"mtgl-std\",\n")
  R_A("\"time\":%le\n", build_time)
  R("},\n")

  free(off); free(ind); free(wgt);
  delete[] src;
  delete[] dst;


  V(Shiloach-Vishkin  Connected components...)

  vertex_property_map<Graph, size_type> componentsMap(g);

  tic();

  size_type count = shiloach_vishkin(g, componentsMap);

  double sv_time = toc();

  R("\"sv\": {\n")
  R("\"name\":\"mtgl-std\",\n")
  R_A("\"time\":%le\n", sv_time)
  R("},\n")

  printf("\tDone %lf\n", sv_time);
  printf("\tComponents %ld\n", count);


  V(BFS...);
  tic();

  vertex_property_map<Graph, size_type> distanceMap(g);

  vertex_iterator verts = vertices(g);
  for(size_type v = 0; v < nv; v++) {
    distanceMap[verts[v]] = 0;
  }

  hop_count_bfs_visitor<Graph, vertex_property_map<Graph, size_type> > hop_count(distanceMap);

  breadth_first_search(g, verts[0], hop_count);

  double sssv_time = toc();

  R("\"sssp\": {\n")
  R("\"name\":\"mtgl-std\",\n")
  R_A("\"time\":%le\n", sssv_time)
  R("},\n")

  printf("\tDone %lf\n", sssv_time);

  V(PageRank...);

  vertex_property_map<Graph, double> ranks(g);
  double epsilon = 1e-8;
  double dampingfactor = 0.85;
  int64_t maxiter = 100;
  tic();

  int iterations = pagerank_omp(g, ranks, epsilon, dampingfactor, maxiter);

  double pr_time = toc();

  R("\"pr\": {\n")
  R("\"name\":\"mtgl-std\",\n")
  R_A("\"time\":%le\n", pr_time)
  R("},\n")

  printf("\tDone %lf\n", pr_time);
  printf("\tIterations %d\n", iterations);

  V(Reading actions...)
  tic();
  
  fp = fopen(argv[2], "r");

  int64_t na;
  int64_t * actions;

  fread(&check, sizeof(uint64_t), 1, fp);

  if(check != endian_check) {
    E(Endianness does not agree.  Order swapping not implemented.);
  }

  fread(&na, sizeof(int64_t), 1, fp);

  actions = (int64_t *)malloc(sizeof(int64_t) * na*2);

  fread(actions, sizeof(int64_t), na*2, fp);

  fclose(fp);

  printf("\t%ld actions read\n", na);

  printf("\tDone %lf\n", toc());

  V(Insert / remove...)
  tic();

  for(uint64_t a = 0; a < na; a++) {
    int64_t i = actions[2*a];
    int64_t j = actions[2*a+1];

    /* is insertion? */
    if(i >= 0) {
      add_edge(verts[i], verts[j], g);
      add_edge(verts[j], verts[i], g);
    } else {
      i = ~i;
      j = ~j;
      //remove_edge(i, j, g);
      //remove_edge(j, i, g);
    }
  }

  double eps = na / toc();

  R("\"update\": {\n")
  R("\"name\":\"mtgl-insertonly\",\n")
  R_A("\"time\":%le\n", eps)
  R("}\n")
  R("},\n")

  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);

  R_A("\"na\":%ld,\n", na)
  R_A("\"mem\":%ld\n", usage.ru_maxrss)
  R("}\n")
  printf("\tDone %lf\n", eps);
  free(actions);
}
