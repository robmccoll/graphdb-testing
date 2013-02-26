#include  <cstdio>

#include  "boost/graph/graph_traits.hpp"
#include  "boost/graph/adjacency_list.hpp"
#include  "boost/graph/undirected_graph.hpp"
#include  "boost/graph/breadth_first_search.hpp"

extern "C" {
#include  "timer.h"
}

#define E_A(X,...) fprintf(stderr, "%s %s %d:\n\t" #X "\n", __FILE__, __func__, __LINE__, __VA_ARGS__); exit(-1);
#define E(X) E_A(X,NULL)
#define V_A(X,...) fprintf(stdout, "%s %s %d:\n\t" #X "\n", __FILE__, __func__, __LINE__, __VA_ARGS__);
#define V(X) V_A(X,NULL)

using namespace boost;

int main(int argc, char *argv[]) {
  if(argc < 3) {
    E_A(Not enough arguments. Usage %s graphfile actionsfile, argv[0]);
  }

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

  V(Creating graph...);

  typedef adjacency_list<vecS, vecS, undirectedS> Graph;

  Graph g(nv);

  V(Loading data into graph...);
  for(uint64_t v = 0; v < nv; v++) {
    for(uint64_t i = off[v]; i < off[v+1]; i++) {
      add_edge(v, ind[i], g);
    }
  }

  free(off); free(ind); free(wgt);

  V(Shiloach-Vishkin  Connected components...)
  int64_t * components = (int64_t *)malloc(sizeof(int64_t) * nv);

  tic();
  for(uint64_t v = 0; v < nv; v++) {
    components[v] = v;
  }

  while(1) {
    uint64_t changed = 0;

    Graph::edge_iterator edgesIt, edgesEnd; tie(edgesIt, edgesEnd) = edges(g);

    for(; edgesIt != edgesEnd; ++edgesIt) {
      if (components[target(*edgesIt,g)] <
	  components[source(*edgesIt,g)]) {
	components[source(*edgesIt,g)] = components[target(*edgesIt,g)];
	changed++;
      }
    }

    if(!changed)
      break;

    for (uint64_t i = 0; i < nv; i++) {
      while (components[i] != components[components[i]])
	components[i] = components[components[i]];
    }
  }

  printf("\tDone %lf\n", toc());
  free(components);

  V(BFS...);
  tic();

  graph_traits<Graph>::vertices_size_type * d = new graph_traits<Graph>::vertices_size_type[nv];
  std::fill_n(d, nv, 0);

  breadth_first_search(g, 0, visitor(make_bfs_visitor(record_distances(d, on_tree_edge()))));

  printf("\tDone %lf\n", toc());

  delete[] d;

  V(PageRank...);

  std::vector<double> tmp_pr(nv);
  std::vector<double> pr(nv);
  double epsilon = 1e-8;
  double dampingfactor = 0.85;
  int64_t maxiter = 100;
  tic();

  std::fill_n(pr.begin(), nv, 1/((double)nv));

  int64_t iter = maxiter;
  double delta = 1;

  while(delta > epsilon && iter > 0) {
    Graph::vertex_iterator vtxIt, vtxEnd; tie(vtxIt, vtxEnd) = vertices(g);
    for(; vtxIt != vtxEnd; ++vtxIt++) {
      tmp_pr[*vtxIt] = 0;

      Graph::out_edge_iterator edgesIt, edgesEnd; 

      tie(edgesIt, edgesEnd) = out_edges(Graph::vertex_descriptor(*vtxIt), g);
      for(; edgesIt != edgesEnd; ++edgesIt) {
	tmp_pr[source(*edgesIt,g)] += (((double)pr[target(*edgesIt, g)]) / 
	  ((double) out_degree(target(*edgesIt, g), g)));
      }
    }

    for(uint64_t v = 0; v < nv; v++) {
      tmp_pr[v] = tmp_pr[v] * dampingfactor + (((double)(1-dampingfactor)) / ((double)nv));
    }

    delta = 0;
    for(uint64_t v = 0; v < nv; v++) {
      double mydelta = tmp_pr[v] - pr[v];

      if(mydelta < 0)
	mydelta = -mydelta;

      delta += mydelta;
      pr[v] = tmp_pr[v];
    }
  }

  printf("\tDone %lf\n", toc());

  pause();

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
      add_edge(i, j, g);
      add_edge(j, i, g);
    } else {
      i = ~i;
      j = ~j;
      remove_edge(i, j, g);
      remove_edge(j, i, g);
    }
  }
  printf("\tDone %lf\n", toc());
  free(actions);
}
