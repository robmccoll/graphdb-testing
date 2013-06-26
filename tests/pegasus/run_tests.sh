#!/bin/bash

VER=1
REDUCERS=40

./run_ccmpt.sh 1024 $REDUCERS graphs/tiny.g.pegasus > tiny.g.cc.$VER.std 2> tiny.g.cc.$VER.err
./run_ccmpt.sh 32768 $REDUCERS graphs/small.g.pegasus > small.g.cc.$VER.std 2> small.g.cc.$VER.err
./run_ccmpt.sh 1048576 $REDUCERS graphs/medium.g.pegasus > medium.g.cc.$VER.std 2> medium.g.cc.$VER.err
./run_ccmpt.sh 16777216 $REDUCERS graphs/large.g.pegasus > large.g.cc.$VER.std 2> large.g.cc.$VER.err

./run_bfs.sh 1024 $REDUCERS graphs/tiny.g.pegasus > tiny.g.sssp.$VER.std 2> tiny.g.sssp.$VER.err
./run_bfs.sh 32768 $REDUCERS graphs/small.g.pegasus > small.g.sssp.$VER.std 2> small.g.sssp.$VER.err
./run_bfs.sh 1048576 $REDUCERS graphs/medium.g.pegasus > medium.g.sssp.$VER.std 2> medium.g.sssp.$VER.err
./run_bfs.sh 16777216 $REDUCERS graphs/large.g.pegasus > large.g.sssp.$VER.std 2> large.g.sssp.$VER.err

./run_pr.sh 1024 $REDUCERS graphs/tiny.g.pegasus nosym > tiny.g.pr.$VER.std 2> tiny.g.pr.$VER.err
./run_pr.sh 32768 $REDUCERS graphs/small.g.pegasus nosym > small.g.pr.$VER.std 2> small.g.pr.$VER.err
./run_pr.sh 1048576 $REDUCERS graphs/medium.g.pegasus nosym > medium.g.pr.$VER.std 2> medium.g.pr.$VER.err
./run_pr.sh 16777216 $REDUCERS graphs/large.g.pegasus nosym > large.g.pr.$VER.std 2> large.g.pr.$VER.err

./run_insert.sh 1024 $REDUCERS graphs/tiny.g.pegasus graphs/tiny.a.pegasus > tiny.g.update.$VER.std 2> tiny.g.update.$VER.err
./run_insert.sh 32768 $REDUCERS graphs/small.g.pegasus graphs/small.a.pegasus > small.g.update.$VER.std 2> small.g.update.$VER.err
./run_insert.sh 1048576 $REDUCERS graphs/medium.g.pegasus graphs/medium.a.pegasus > medium.g.update.$VER.std 2> medium.g.update.$VER.err
./run_insert.sh 16777216 $REDUCERS graphs/large.g.pegasus graphs/large.a.pegasus > large.g.update.$VER.std 2> large.g.update.$VER.err

