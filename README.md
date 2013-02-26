Graph Package Testing
=====================

The goal of this project is to benchmark various graph databases, engines, datastructures, 
and data stores. The packages below will be measured in terms of performance on these graph 
algorithms which are all fairly basic and cover a reasonable range of styles (data structure
modification, breadth first traversal, edge-parallel, and bulk-synchronous parallel).

- Insertion / deletion / update 
- Single-source shortest paths (output must be unweighted distance)
- Shiloach-Vishkin style connected components
- PageRank

Additionally, information about their licensing, costs, capabilities, distribution patterns,
etc. will be gathered but may or may not be available here.

| Package Name  | Type       |S-V Components| SSSP-BFT(BFS)| PageRank     |Insert/Remove |
|---------------|:----------:|:------------:|:------------:|:------------:|:------------:|
| STINGER       | Library    | Implemented  | Implemented  | Implemented  | Implemented  |
| MySQL         | SQLDB      |              |              |              |              |
| SQLite        | SQLDB      | Implemented  | Implemented  | Implemented  | Implemented  |
| Oracle        | SQLDB      |              |              |              |              |
| SQL Server    | SQLDB      |              |              |              |              |
| Neo4j         | GDB/NoSQL  |              |              |              |              |
| OrientDB      | GDB/NoSQL  |              |              |              |              |
| InfoGrid      | GDB        |              |              |              |              |
| Titan         | GDB        |              |              |              |              |
| FlockDB       | GDB        |              |              |              |              |
| ArangoDB      | GDB/KV/DOC |              |              |              |              |
| InfiniteGraph | GDB        |              |              |              |              |
| AllegroGraph  | GDB        |              |              |              |              |
| DEX           | GDB        |              |              |              |              |
| GraphBase     | GDB        |              |              |              |              |
| HyperGraphDB  | HyperGDB   |              |              |              |              |
| Bagel         | BSP        |              |              |              |              |
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
                                             
