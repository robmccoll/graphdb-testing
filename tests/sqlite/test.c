#include    "sqlite3.h"
#include    "timer.h"

#include  <stdio.h>
#include  <stdlib.h>
#include  <stdint.h>
#include <sys/time.h>
#include <sys/resource.h>

#define E_A(X,...) fprintf(stderr, "%s %s %d:\n\t" #X "\n", __FILE__, __func__, __LINE__, __VA_ARGS__); exit(-1);
#define E(X) E_A(X,NULL)
#define V_A(X,...) fprintf(stdout, "%s %s %d:\n\t" #X "\n", __FILE__, __func__, __LINE__, __VA_ARGS__);
#define V(X) V_A(X,NULL)
#define R_A(X,...) fprintf(stdout, "RSLT: " X, __VA_ARGS__);
#define R(X) R_A(X,NULL)

#define DB_OR_DIE(X) \
  if(SQLITE_OK != sqlite3_exec(db, \
    X, NULL, 0, &zErrMsg)) { \
    E_A(SQL command failed: %s\nError: %s, X, zErrMsg); \
  }

static int
print_result(void *NotUsed, int argc, char ** argv, char **azColName) {
  for(int i = 0; i < argc; i++){
    printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
  }
  printf("\n");
  return 0;
}

static int
get_count(void * count, int argc, char ** argv, char **azColName) {
  *((int64_t *)count) = atol(argv[0]);
  return 0;
}

static int
get_double(void * rtn, int argc, char ** argv, char **azColName) {
  *((double *)rtn) = atof(argv[0]);
  return 0;
}

int main(int argc, char *argv[]) {
  if(argc < 3) {
    E_A(Not enough arguments. Usage %s graphfile actionsfile, argv[0]);
  }

  R("{\n")
  R("\"type\":\"sqlite\",\n")

  int64_t iter = 0;

  V(Loading graph...);

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

  off = malloc(sizeof(int64_t) * nv+1);
  ind = malloc(sizeof(int64_t) * ne);
  wgt = malloc(sizeof(int64_t) * ne);

  fread(off, sizeof(int64_t), nv+1, fp);
  fread(ind, sizeof(int64_t), ne, fp);
  fread(wgt, sizeof(int64_t), ne, fp);

  fclose(fp);

  R_A("\"nv\":%ld,\n", nv)
  R_A("\"ne\":%ld,\n", ne)
  R("\"results\": {\n")

  V(Creating db and edges table...);

  sqlite3 *db;
  char *zErrMsg = 0;
  int rc;

  if(getenv("DB_FILE")) {
    if(sqlite3_open(getenv("DB_FILE"), &db)) {
      sqlite3_close(db);
      E(failed to open db);
    }
  } else {
    if(sqlite3_open(":memory:", &db)) {
      sqlite3_close(db);
      E(failed to open db);
    }
  }

  if(SQLITE_OK != sqlite3_exec(db, 
    "DROP TABLE IF EXISTS edges", NULL, 0, &zErrMsg)) {
    E(Creating edge table failed);
  }

  if(SQLITE_OK != sqlite3_exec(db, 
    "CREATE TABLE edges (src BIG INT NOT NULL, dst BIG INT NOT NULL, wgt BIT INT NOT NULL)", NULL, 0, &zErrMsg)) {
    E(Creating edge table failed);
  }
  DB_OR_DIE("CREATE UNIQUE INDEX IF NOT EXISTS edgepairs ON edges (src, dst)");
  DB_OR_DIE("CREATE INDEX IF NOT EXISTS edgesrcs ON edges (src)");
  DB_OR_DIE("CREATE INDEX IF NOT EXISTS edgedsts ON edges (dst)");

  V(Loading data into edges table...);
  char sqlcmd[1024];

  tic();
  for(uint64_t v = 0; v < nv; v++) {
    for(uint64_t i = off[v]; i < off[v+1]; i++) {
      sprintf(sqlcmd, "INSERT OR IGNORE INTO edges (src, dst, wgt) VALUES (%ld, %ld, %ld)", v, ind[i], wgt[i]);
      DB_OR_DIE(sqlcmd);
    }
  }
  double build_time = toc();
  R("\"build\": {\n")
  R("\"name\":\"sqlite-std\",\n")
  R_A("\"time\":%le\n", build_time)
  R("},\n")

  free(off); free(ind); free(wgt);

  V(Setting up connected components...);
  tic();

  DB_OR_DIE("DROP TABLE IF EXISTS components");
  DB_OR_DIE("CREATE TABLE components (vtx BIGINT UNIQUE, label BIGINT)");
  DB_OR_DIE("CREATE UNIQUE INDEX IF NOT EXISTS pairs ON components (vtx, label)");
  DB_OR_DIE("CREATE UNIQUE INDEX IF NOT EXISTS src ON components (vtx)");
  DB_OR_DIE("CREATE INDEX IF NOT EXISTS labels ON components (label)");

  DB_OR_DIE("DROP TABLE IF EXISTS components_new");
  DB_OR_DIE("CREATE TABLE components_new (vtx BIGINT UNIQUE, label BIGINT)");
  DB_OR_DIE("CREATE UNIQUE INDEX IF NOT EXISTS pairs ON components_new (vtx, label)");
  DB_OR_DIE("CREATE UNIQUE INDEX IF NOT EXISTS src ON components_new (vtx)");
  DB_OR_DIE("CREATE INDEX IF NOT EXISTS labels ON components_new (label)");
  
  printf("\tDone %lf\n", toc());

  V(Performing connected components...);
  tic();

  for(uint64_t v = 0; v < nv; v++) {
      sprintf(sqlcmd, "INSERT INTO components (vtx, label) VALUES (%ld, %ld)", v, v);
      DB_OR_DIE(sqlcmd);
      sprintf(sqlcmd, "INSERT INTO components_new (vtx, label) VALUES (%ld, %ld)", v, v);
      DB_OR_DIE(sqlcmd);
  }

  uint64_t old_count = nv;
  iter = 0;
  while(1) {
    DB_OR_DIE("UPDATE components_new SET label = ( "
        "SELECT MIN(CASE WHEN src.label > dst.label THEN dst.label ELSE src.label END) "
        "FROM edges "
        "LEFT JOIN components AS src ON edges.src = src.vtx "
        "LEFT JOIN components AS dst ON edges.dst = dst.vtx "
        "WHERE components_new.vtx = src.vtx "
        "GROUP BY edges.src "
        ")");

    DB_OR_DIE("DELETE FROM components");
    DB_OR_DIE("INSERT INTO components SELECT components_new.vtx, IFNULL(components_new.label, components_new.vtx)  FROM components_new");

    int64_t new_count = 0;
    sqlite3_exec(db, "SELECT COUNT(DISTINCT label) FROM components", get_count, &new_count, &zErrMsg);

    printf("\tIteration %ld: Count is %ld\n", iter, new_count);
    if(old_count == new_count)
      break;
    else
      old_count = new_count;
  }

  double sv_time = toc();

  R("\"sv\": {\n")
  R("\"name\":\"sqlite-std\",\n")
  R_A("\"time\":%le\n", sv_time)
  R("},\n")

  printf("\tDone %lf\n", sv_time);

  V(Setting up BFS...);
  tic();

  DB_OR_DIE("DROP TABLE IF EXISTS distance");
  DB_OR_DIE("CREATE TABLE distance (vtx BIGINT NOT NULL, dist BIGINT NOT NULL)");
  DB_OR_DIE("CREATE UNIQUE INDEX IF NOT EXISTS nv ON distance (vtx)");

  printf("\tDone %lf\n", toc());

  V(Performing BFS...);
  tic();

  uint64_t dist = 0;
  uint64_t start = nv / 4;
  sprintf(sqlcmd, "INSERT INTO distance (vtx, dist) VALUES (%ld, %ld)", start, dist);
  DB_OR_DIE(sqlcmd);

  while (1) {
    uint64_t newdist = dist + 1;
    sprintf(sqlcmd, 
      "INSERT OR IGNORE INTO distance (vtx, dist) "
      "SELECT DISTINCT edges.dst, %ld "
      "FROM edges "
      "LEFT JOIN distance "
      "ON edges.src = distance.vtx "
      "WHERE distance.dist = %ld", newdist, dist);
    DB_OR_DIE(sqlcmd);

    int rows_affected = sqlite3_changes(db);
    printf("\tAffected %d rows\n", rows_affected);

    if(rows_affected < 1)
      break;

    dist++;
  }

  double sssv_time = toc();

  R("\"sssp\": {\n")
  R("\"name\":\"sqlite-std\",\n")
  R_A("\"time\":%le\n", sssv_time)
  R("},\n")

  printf("\tDone %lf\n", sssv_time);

  V(Setting up PageRank...);
  tic();

  DB_OR_DIE("DROP TABLE IF EXISTS outdegree");
  DB_OR_DIE("CREATE TABLE outdegree (vtx BIGINT NOT NULL, degree BIGINT NOT NULL) ");
  DB_OR_DIE("CREATE INDEX IF NOT EXISTS deg ON outdegree (vtx, degree)");

  DB_OR_DIE("DROP TABLE IF EXISTS pagerank");
  DB_OR_DIE("CREATE TABLE pagerank (vtx BIGINT NOT NULL, pagerank DOUBLE NOT NULL) ");
  DB_OR_DIE("CREATE UNIQUE INDEX IF NOT EXISTS nv ON pagerank (vtx)");

  DB_OR_DIE("DROP TABLE IF EXISTS pagerank_new");
  DB_OR_DIE("CREATE TABLE pagerank_new (vtx BIGINT NOT NULL, pagerank DOUBLE NOT NULL) ");
  DB_OR_DIE("CREATE UNIQUE INDEX IF NOT EXISTS nv ON pagerank_new (vtx)");

  /* Give this for free - we could count this as part of datastructure init */
  DB_OR_DIE("INSERT INTO outdegree (vtx, degree) SELECT src, COUNT(src) FROM edges GROUP BY src");

  printf("\tDone %lf\n", toc());

  V(Performing PageRank...);
  tic();

  double startPR = 1.0 / ((double)nv);
  sprintf(sqlcmd, "INSERT INTO pagerank (vtx, pagerank) SELECT DISTINCT src, %lf FROM edges", startPR);
  DB_OR_DIE(sqlcmd);

  double delta = 1.0;
  double epsilon = 1e-8;
  double dampingfactor = 0.85;
  double damping = (1.0 - dampingfactor) / ((double)nv);

  uint64_t maxiter = 100;
  iter = maxiter;

  while (delta > epsilon && iter > 0) {
    DB_OR_DIE("INSERT INTO pagerank_new (vtx, pagerank) "
	      "SELECT edges.src, SUM(pagerank/degree) AS newPR "
	      "FROM edges "
	      "LEFT JOIN pagerank ON edges.dst = pagerank.vtx "
	      "LEFT JOIN outdegree ON edges.dst = outdegree.vtx "
	      "GROUP BY edges.src");

    sprintf(sqlcmd, "UPDATE pagerank_new SET pagerank = %lf * pagerank + %lf", dampingfactor, damping);
    DB_OR_DIE(sqlcmd);

    sqlite3_exec(db,"SELECT SUM(ABS(pagerank.pagerank-pagerank_new.pagerank)) AS delta "
	      "FROM pagerank "
	      "LEFT JOIN pagerank_new ON pagerank.vtx = pagerank_new.vtx;", get_double, &delta, &zErrMsg);

    DB_OR_DIE("DELETE FROM pagerank");
    DB_OR_DIE("INSERT INTO pagerank SELECT * FROM pagerank_new");
    DB_OR_DIE("DELETE FROM pagerank_new");

    printf("\tdelta: %lf\n", delta);
    printf("\titer: %ld\n", maxiter - iter + 1);
    iter--;
  }

  double pr_time = toc();

  R("\"pr\": {\n")
  R("\"name\":\"sqlite-std\",\n")
  R_A("\"time\":%le\n", pr_time)
  R("},\n")

  printf("\tDone %lf\n", pr_time);

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

  actions = malloc(sizeof(int64_t) * na*2);

  fread(actions, sizeof(int64_t), na*2, fp);

  fclose(fp);

  printf("\t%ld actions read\n", na);

  printf("\tDone %lf\n", toc());

  V(Insert remove test...)
  tic();

  for(uint64_t a = 0; a < na; a++) {
    int64_t i = actions[2*a];
    int64_t j = actions[2*a+1];

    /* is insertion? */
    if(i >= 0) {
      sprintf(sqlcmd, "INSERT OR IGNORE INTO edges (src, dst, wgt) VALUES (%ld, %ld, 1) ", i, j);
      DB_OR_DIE(sqlcmd);
      sprintf(sqlcmd, "INSERT OR IGNORE INTO edges (src, dst, wgt) VALUES (%ld, %ld, 1) ", j, i);
      DB_OR_DIE(sqlcmd);
    } else {
      i = ~i;
      j = ~j;

      sprintf(sqlcmd, "DELETE FROM edges WHERE src = %ld AND dst = %ld", i, j);
      DB_OR_DIE(sqlcmd);
      sprintf(sqlcmd, "DELETE FROM edges WHERE src = %ld AND dst = %ld", j, i);
      DB_OR_DIE(sqlcmd);
    }
  }

  double eps = na / toc();

  R("\"update\": {\n")
  R("\"name\":\"sqlite-std\",\n")
  R_A("\"time\":%le\n", eps)
  R("}\n")
  R("},\n")

  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);

  R_A("\"na\":%ld,\n", na)
  R_A("\"mem\":%ld\n", usage.ru_maxrss)
  R("}\n")

  printf("\tDone %lf\n", eps);

  V(Closing db...);
  sqlite3_close(db);
  V(Exiting);

  return 0;
}
