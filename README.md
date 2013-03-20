Graph Package Testing
=====================

The goal of this project is to benchmark various graph databases, engines, datastructures, 
and data stores. The packages below will be measured in terms of performance on these graph 
algorithms which are all fairly basic and cover a reasonable range of styles (data structure
modification, breadth first traversal, edge-parallel, and bulk-synchronous parallel).

- Insertion / deletion / update 
  - Implemented according to the standard semantics of the data structure
- Single-source shortest paths (output must be unweighted distance)
  - Must be implemented using a breadth-first traversal
- Shiloach-Vishkin style connected components
  - Must be implemented using a (parallel where possible) for all edges loop
- PageRank
  - Must be implemented using (parallel where possible) for all vertices loops

Additionally, information about their licensing, costs, capabilities, distribution patterns,
etc. will be gathered but may or may not be available here.

Where possible, libraries will be pulled in directly from their git repositories 
as submodules.  In cases where submodules cannot be used, the source / binary necessary
for the library will be added to the repository if the license permits.  Other cases
should provide README files or something similar where the files should be placed and 
how to obtain them.  For the most part, projects which are not available as free and/or
open source software will only be included or tested if they have an evaluation license
of some sort.  In those cases, only the evaluation version is used.  For example, DEX
has a limitation on the evaluation version that prevents more than 1M total vertices and
edges from being added to the graph which is the version used here.


| Package Name  | Type       |S-V Components| SSSP-BFT(BFS)| PageRank     |Insert/Remove |
|---------------|:----------:|:------------:|:------------:|:------------:|:------------:|
| STINGER       | Library    | Implemented  | Implemented  | Implemented  | Implemented  |
| MySQL         | SQLDB      |              |              |              |              |
| SQLite        | SQLDB      | Implemented  | Implemented  | Implemented  | Implemented  |
| Oracle        | SQLDB      | Not Planned  | Not Planned  | Not Planned  | Not Planned  |
| SQL Server    | SQLDB      | Not Planned  | Not Planned  | Not Planned  | Not Planned  |
| Neo4j         | GDB/NoSQL  | Implemented  | Implemented  | Implemented  | Implemented  |
| OrientDB	| GDB/NoSQL  | Implemented  | Implemented  | Implemented  | Implemented  |
| InfoGrid      | GDB        |              |              |              |              |
| Titan         | GDB        | Implemented  | Implemented  | Implemented  | Implemented  |
| FlockDB       | GDB        |              |              |              |              |
| ArangoDB      | GDB/KV/DOC | Implemented  | Implemented  | Implemented  |              |
| InfiniteGraph | GDB        |              |              |              |              |
| AllegroGraph  | GDB        |              |              |              |              |
| DEX           | GDB        | Implemented  | Implemented  | Implemented  | Implemented  |
| GraphBase     | GDB        |              |              |              |              |
| HyperGraphDB  | HyperGDB   |              |              |              |              |
| Bagel         | BSP        | Implemented  | Implemented  | Implemented  | Implemented  |
| Hama          | BSP        |              |              |              |              |
| Giraph        | BSP        |              |              |              |              |
| PEGASUS       | Hadoop     |              |              |              |              |
| Faunus        | Hadoop     |              |              |              |              |
| NetworkX      | Library    | Implemented  | Implemented  | Implemented  | Implemented  |
| Gephi         | Toolkit    |              |              |              |              |
| MTGL          | Library    | Implemented  | Implemented  | Implemented  | Implemented  |
| BGL           | Library    | Implemented  | Implemented  | Implemented  |              |
| GraphStream   | Library    |              |              |              |              |
| uRika         | Appliance  | N/A          | N/A          | N/A          | N/A          |
                                             
Graphs will be generated using an implementation of the [R-MAT](http://repository.cmu.edu/compsci/541/)
synthetic graph generator which is designed to generate graphs that emulate the properties of real
social networks at a large scale (small-world phenomena, power-law degree distribution, etc).
Graphs are generated using parameters A = 0.55, B = 0.1, C = 0.1, D = 0.25.  The size of the
graph is determined by SCALE and EDGE FACTOR, where the number of vertices = 2^SCALE and the number of
edges is the number of vertices multiplied by the edge factor.  Inserted, updated, and deleted edges
are generated following this same procedure with P(delete) = 0.0625 that instead of generating an 
insertion, a previous inserted edge will be selected for deletion. Graph sizes used for testing
are listed below:

| Name    | SCALE | EDGE FACTOR | Vertices | Edges |
|---------|-------|-------------|----------|-------|
| Tiny    | 10    | 8           | 1K       | 8K    |
| Small   | 15    | 8           | 32K      | 256K  |
| Medium  | 20    | 8           | 1M       | 8M    |
| Large   | 24    | 8           | 16M      | 128M  |
