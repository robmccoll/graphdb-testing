#include "gdb/Dex.h"
#include "gdb/Database.h"
#include "gdb/Session.h"
#include "gdb/Graph.h"
#include "gdb/Objects.h"
#include "gdb/ObjectsIterator.h"
#include "algorithms/TraversalBFS.h"

#include  <stdio.h>
#include  <algorithm>
#include <sys/time.h>
#include <sys/resource.h>
#include  <map>

extern "C" {
#include  "timer.h"
}

#define E_A(X,...) fprintf(stderr, "%s %s %d:\n\t" #X "\n", __FILE__, __func__, __LINE__, __VA_ARGS__); exit(-1);
#define E(X) E_A(X,NULL)
#define V_A(X,...) fprintf(stdout, "%s %s %d:\n\t" #X "\n", __FILE__, __func__, __LINE__, __VA_ARGS__);
#define V(X) V_A(X,NULL)
#define R_A(X,...) fprintf(stdout, "RSLT: " X, __VA_ARGS__);
#define R(X) R_A(X,NULL)

using namespace dex::gdb;
using namespace dex::algorithms;

/* to fix ambiguity issues - its this or be totally verbose on 
 * namespaces. c++ is fun */
typedef long int int64;
typedef long unsigned int uint64;

int main(int argc, char *argv[])
{
  if(argc < 3) {
    E_A(Not enough arguments. Usage %s graphfile actionsfile, argv[0]);
  }

  R("{\n")
  R("\"type\":\"dex\",\n")

  FILE * fp = fopen(argv[1], "r");

  int64 nv;
  int64 ne;
  int64 * off, * ind, * wgt;

  const int64 endian_check = 0x1234ABCDul;
  int64 check;

  fread(&check, sizeof(int64), 1, fp);

  if(check != endian_check) {
    E(Endianness does not agree.  Order swapping not implemented.);
  }

  fread(&nv, sizeof(int64), 1, fp);
  fread(&ne, sizeof(int64), 1, fp);

  off = (int64 *)malloc(sizeof(int64) * nv+1);
  ind = (int64 *)malloc(sizeof(int64) * ne);
  wgt = (int64 *)malloc(sizeof(int64) * ne);

  fread(off, sizeof(int64), nv+1, fp);
  fread(ind, sizeof(int64), ne, fp);
  fread(wgt, sizeof(int64), ne, fp);

  fclose(fp);

  R_A("\"nv\":%ld,\n", nv)
  R_A("\"ne\":%ld,\n", ne)
  R("\"results\": {\n")

  V(Creating graph...);

  oid_t * vertices = new oid_t[nv];
  std::fill_n(vertices, nv, -1);

  tic();
  DexConfig cfg;
  Dex * dex = new Dex(cfg);
  Database * db = dex->Create(L"./graphdb.dex", L"GraphDB");
  Session * sess = db->NewSession();
  Graph * graph = sess->GetGraph();

  type_t vtxType = graph->NewNodeType(L"Vertex");

  type_t edgeType = graph->NewEdgeType(L"Edge", false, false);
  attr_t edgeWeightType = graph->NewAttribute(edgeType, L"weight", Long, Basic);

  Value * value = new Value();

  for(int64 v = 0; v < nv; v++) {
    if(vertices[v] == -1) {
      vertices[v] = graph->NewNode(vtxType);
    }
    for(int64 i = off[v]; i < off[v+1]; i++) {
      int64 u = ind[i];
      if(vertices[u] == -1) {
	vertices[u] = graph->NewNode(vtxType);
      }
      oid_t edge = graph->NewEdge(edgeType, vertices[u], vertices[v]);
      graph->SetAttribute(edge, edgeWeightType, value->SetLong(wgt[i]));
    }
  }

  double build_time = toc();
  R("\"build\": {\n")
  R("\"name\":\"dex-std\",\n")
  R_A("\"time\":%le\n", build_time)
  R("},\n")

#if 1
  V(Shiloach-Vishkin  Connected components...)
  std::map<oid_t, int64> components;

  tic();
  {
    Objects * vtxObjects = graph->Select(vtxType);
    ObjectsIterator * vtx = vtxObjects->Iterator();
    while(vtx->HasNext()) {
      oid_t cur = vtx->Next();
      components[cur] = cur;
    }
    delete vtx;
    delete vtxObjects;
  }

  Objects * edgeObjects = graph->Select(edgeType);
  while(1) {
    uint64 changed = 0;

    ObjectsIterator * edges = edgeObjects->Iterator();

    while(edges->HasNext()) {
      oid_t cur = edges->Next();
      EdgeData * data = graph->GetEdgeData(cur);
      if (components[data->GetTail()] <
	  components[data->GetHead()]) {
	components[data->GetHead()] = components[data->GetTail()];
	changed++;
      }
      delete data;
    }

    delete edges;

    if(!changed)
      break;

    for (uint64 i = 0; i < nv; i++) {
      while (components[i] != components[components[i]])
	components[i] = components[components[i]];
    }
  }
  delete edgeObjects;

  double sv_time = toc();

  R("\"sv\": {\n")
  R("\"name\":\"dex-std\",\n")
  R_A("\"time\":%le\n", sv_time)
  R("},\n")

  printf("\tDone %lf\n", sv_time);

  V(BFS...);
  tic();

  std::map<oid_t, int64> distance;
  distance[vertices[0]] = 0;
  
  {
    TraversalBFS bfs(*sess, vertices[0]);
    bfs.AddAllEdgeTypes(Outgoing);
    bfs.AddAllNodeTypes();
    while(bfs.HasNext()) {
      distance[bfs.Next()] = bfs.GetCurrentDepth();
    }
  }


  double sssv_time = toc();

  R("\"sssp\": {\n")
  R("\"name\":\"dex-std\",\n")
  R_A("\"time\":%le\n", sssv_time)
  R("},\n")

  printf("\tDone %lf\n", sssv_time);

  V(PageRank...);

  std::map<oid_t,double> tmp_pr;
  std::map<oid_t,double> pr;
  double epsilon = 1e-8;
  double dampingfactor = 0.85;
  int64 maxiter = 100;
  tic();

  Objects * vtxObjects = graph->Select(vtxType);
  {
    ObjectsIterator * vtx = vtxObjects->Iterator();
    while(vtx->HasNext()) {
      oid_t cur = vtx->Next();
      pr[cur] = 1/((double)nv);
    }
    delete vtx;
  }

  int64 iter = maxiter;
  double delta = 1;

  while(delta > epsilon && iter > 0) {
    {
      ObjectsIterator * vtx = vtxObjects->Iterator();
      while(vtx->HasNext()) {
	oid_t cur = vtx->Next();
	tmp_pr[cur] = 0;

	Objects * neighborObjects = graph->Neighbors(cur, edgeType, Outgoing);
	ObjectsIterator * neighbor = neighborObjects->Iterator();
	while(neighbor->HasNext()) {
	  oid_t neigh = neighbor->Next();
	  tmp_pr[cur] += (pr[neigh] / ((double) graph->Degree(neigh, edgeType, Outgoing)));
	}
	delete neighbor;
	delete neighborObjects;
      }
      delete vtx;
    }

    {
      ObjectsIterator * vtx = vtxObjects->Iterator();
      while(vtx->HasNext()) {
	oid_t cur = vtx->Next();
	tmp_pr[cur] = tmp_pr[cur] * dampingfactor + (((double)(1-dampingfactor))/((double)nv));
      }
      delete vtx;
    }

    delta = 0;
    {
      ObjectsIterator * vtx = vtxObjects->Iterator();
      while(vtx->HasNext()) {
	oid_t cur = vtx->Next();
	double mydelta = tmp_pr[cur] - pr[cur];

	if(mydelta < 0)
	  mydelta = -mydelta;

	delta += mydelta;
	pr[cur] = tmp_pr[cur];
      }
      delete vtx;
    }

    iter--;
  }

  delete vtxObjects;

  double pr_time = toc();

  R("\"pr\": {\n")
  R("\"name\":\"dex-std\",\n")
  R_A("\"time\":%le\n", pr_time)
  R("},\n")

  printf("\tDone %lf\n", pr_time);
#endif

  V(Reading actions...)
  tic();
  
  fp = fopen(argv[2], "r");

  int64 na;
  int64 * actions;

  fread(&check, sizeof(uint64), 1, fp);

  if(check != endian_check) {
    E(Endianness does not agree.  Order swapping not implemented.);
  }

  fread(&na, sizeof(int64), 1, fp);

  actions = (int64 *)malloc(sizeof(int64) * na*2);

  fread(actions, sizeof(int64), na*2, fp);

  fclose(fp);

  printf("\t%ld actions read\n", na);

  printf("\tDone %lf\n", toc());

  V(Insert / remove...)
  tic();

  if(na > 100000) {
    V(Due to object limit lowering insert / delete to 100k);
    na = 100000;
  }

  for(uint64 a = 0; a < na; a++) {
    int64 i = actions[2*a];
    int64 j = actions[2*a+1];

    /* is insertion? */
    if(i >= 0) {
      if(vertices[i] == -1) {
	vertices[i] = graph->NewNode(vtxType);
      }
      if(vertices[j] == -1) {
	vertices[j] = graph->NewNode(vtxType);
      }
      oid_t edge = graph->FindEdge(edgeType, vertices[i], vertices[j]);
      if(Objects::InvalidOID !=  edge) {
	graph->GetAttribute(edge, edgeWeightType, *value);
	value->SetLong(value->GetLong()+1);
	graph->SetAttribute(edge, edgeWeightType, *value);
      } else {
	edge = graph->NewEdge(edgeType, vertices[i], vertices[j]);
	graph->SetAttribute(edge, edgeWeightType, value->SetLong(1));
      }
      edge = graph->FindEdge(edgeType, vertices[j], vertices[i]);
      if(Objects::InvalidOID !=  edge) {
	graph->GetAttribute(edge, edgeWeightType, *value);
	value->SetLong(value->GetLong()+1);
	graph->SetAttribute(edge, edgeWeightType, *value);
      } else {
	edge = graph->NewEdge(edgeType, vertices[j], vertices[i]);
	graph->SetAttribute(edge, edgeWeightType, value->SetLong(1));
      }
    } else {
      i = ~i;
      j = ~j;
      if(vertices[i] != -1 && vertices[j] != -1) {
	oid_t edge = graph->FindEdge(edgeType, vertices[i], vertices[j]);
	if(Objects::InvalidOID !=  edge) {
	  graph->Drop(edge);
	}
      }
    }
  }

  double eps = na / toc();

  R("\"update\": {\n")
  R("\"name\":\"dex-std\",\n")
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

  delete value;


  delete sess;
  delete db;
  delete dex;

  return 0;
}

