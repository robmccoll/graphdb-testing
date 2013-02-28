ARANGODB
========

For ArangoDB, I'm doing something a bit different.  The static graph loading and input / output
are implemented in php and the connected components, BFS, and PageRank are implemented in 
JavaScript with support from AQL.  php is required for the file binary file parsing. 
JavaScript is used since it is the default interface to ArangoDB and the interface is
much simpler than its php counterpart.

Runing the connected components, BFS, and PageRank should be as easy as (or something similar):
arangosh --server.disable-authentication true --javascript.execute test.js 
