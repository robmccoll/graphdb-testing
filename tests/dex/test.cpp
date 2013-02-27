#include "gdb/Dex.h"
#include "gdb/Database.h"
#include "gdb/Session.h"
#include "gdb/Graph.h"
#include "gdb/Objects.h"
#include "gdb/ObjectsIterator.h"
#include "algorithms/TraversalBFS.h"

#include  <stdio.h>
#include  <algorithm>
#include  <map>

extern "C" {
#include  "timer.h"
}

#define E_A(X,...) fprintf(stderr, "%s %s %d:\n\t" #X "\n", __FILE__, __func__, __LINE__, __VA_ARGS__); exit(-1);
#define E(X) E_A(X,NULL)
#define V_A(X,...) fprintf(stdout, "%s %s %d:\n\t" #X "\n", __FILE__, __func__, __LINE__, __VA_ARGS__);
#define V(X) V_A(X,NULL)

using namespace dex::gdb;
using namespace dex::algorithms;

typedef long int int64;
typedef long unsigned int uint64;

int main(int argc, char *argv[])
{
  if(argc < 3) {
    E_A(Not enough arguments. Usage %s graphfile actionsfile, argv[0]);
  }

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

  V(Creating graph...);

  oid_t * vertices = new oid_t[nv];
  std::fill_n(vertices, nv, -1);

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

  delete value;

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

  while(1) {
    uint64 changed = 0;

    Objects * edgeObjects = graph->Select(edgeType);
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
    delete edgeObjects;

    if(!changed)
      break;

    for (uint64 i = 0; i < nv; i++) {
      while (components[i] != components[components[i]])
	components[i] = components[components[i]];
    }
  }

  printf("\tDone %lf\n", toc());

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

  printf("\tDone %lf\n", toc());

  delete sess;
  delete db;
  delete dex;

  return 0;
}

