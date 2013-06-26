#!/bin/bash

VER=3
WORKERS=16

GRAPH_DIR=graphs
GRAPH=tiny

for GRAPH in tiny small medium large
  do

  time ./run_sssp.sh $GRAPH.g $WORKERS > $GRAPH.g.sssp.$VER.std 2> $GRAPH.g.sssp.$VER.err
  time ./run_cc.sh $GRAPH.g $WORKERS > $GRAPH.g.cc.$VER.std 2> $GRAPH.g.cc.$VER.err
  time ./run_pr.sh $GRAPH.g $WORKERS > $GRAPH.g.pr.$VER.std 2> $GRAPH.g.pr.$VER.err
  time ./run_insertremove.sh $GRAPH.g /home/dediger/$GRAPH.a.pegasus $WORKERS > $GRAPH.g.insert.$VER.std 2> $GRAPH.g.insert.$VER.err

done
