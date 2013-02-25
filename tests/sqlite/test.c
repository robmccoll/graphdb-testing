#include    "sqlite3.h"

#include  <stdio.h>
#include  <stdlib.h>
#include  <stdint.h>

#define E_A(X,...) fprintf(stderr, "%s %s %d:\n\t" #X "\n", __FILE__, __func__, __LINE__, __VA_ARGS__); exit(-1);
#define E(X) E_A(X,NULL)
#define V_A(X,...) fprintf(stdout, "%s %s %d:\n\t" #X "\n", __FILE__, __func__, __LINE__, __VA_ARGS__);
#define V(X) V_A(X,NULL)

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

int main(int argc, char *argv[]) {
  if(argc < 3) {
    E_A(Not enough arguments. Usage %s graphfile actionsfile, argv[0]);
  }

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
    "CREATE TABLE edges (src BIG INT NOT NULL, dst BIG INT NOT NULL)", NULL, 0, &zErrMsg)) {
    E(Creating edge table failed);
  }

  V(Loading data into edges table...);
  char sqlcmd[1024];

  for(uint64_t v = 0; v < nv; v++) {
    for(uint64_t i = off[v]; i < off[v+1]; i++) {
      sprintf(sqlcmd, "INSERT INTO edges (src, dst) VALUES (%ld, %ld)", v, ind[i]);
      DB_OR_DIE(sqlcmd);
    }
  }

  free(off); free(ind); free(wgt);

  V(Setting up connected components...);

  DB_OR_DIE("DROP TABLE IF EXISTS components");
  DB_OR_DIE("CREATE TABLE components (vtx BIGINT NOT NULL, label BIGINT NOT NULL)");
  DB_OR_DIE("CREATE UNIQUE INDEX IF NOT EXISTS pairs ON components (vtx)");

  DB_OR_DIE("DROP TABLE IF EXISTS components_new");
  DB_OR_DIE("CREATE TABLE components_new (vtx BIGINT NOT NULL, label BIGINT NOT NULL)");
  DB_OR_DIE("CREATE UNIQUE INDEX IF NOT EXISTS pairs ON components_new (vtx)");

  for(uint64_t v = 0; v < nv; v++) {
      sprintf(sqlcmd, "INSERT INTO components (vtx, label) VALUES (%ld, %ld)", v, v);
      DB_OR_DIE(sqlcmd);
      sprintf(sqlcmd, "INSERT INTO components_new (vtx, label) VALUES (%ld, %ld)", v, v);
      DB_OR_DIE(sqlcmd);
  }

  while(1) {
    DB_OR_DIE("UPDATE components_new SET label = IFNULL(( "
        "SELECT MIN(CASE WHEN src.label > dst.label THEN dst.label ELSE src.label END) "
        "FROM edges "
        "LEFT JOIN components AS src ON edges.src = src.vtx "
        "LEFT JOIN components AS dst ON edges.dst = dst.vtx "
        "WHERE components_new.vtx = src.vtx AND src.label != dst.label "
        "GROUP BY edges.src "
        "), vtx)");

    int rows_affected = sqlite3_changes(db);
    printf("\tAffected %d rows\n", rows_affected);

    if(rows_affected < 1)
      break;

    DB_OR_DIE("DELETE FROM components");
    DB_OR_DIE("INSERT INTO components SELECT * FROM components_new");
  }

  if(SQLITE_OK != sqlite3_exec(db, 
    "SELECT COUNT(DISTINCT `label`) AS numComp FROM components", print_result, 0, &zErrMsg)) {
    E(Creating edge table failed);
  }

  V(Setting up BFS...);

  DB_OR_DIE("DROP TABLE IF EXISTS distance");
  DB_OR_DIE("CREATE TABLE distance (vtx BIGINT NOT NULL, dist BIGINT NOT NULL)");
  DB_OR_DIE("CREATE UNIQUE INDEX IF NOT EXISTS nv ON distance (vtx)");

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

  V(Reading actions...)
  
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

  V(Insert remove test...)

  //TODO

  V(Closing db...);
  sqlite3_close(db);
  V(Exiting);
  return 0;
}
