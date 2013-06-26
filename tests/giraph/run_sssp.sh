#!/bin/bash

WORKERS=$2

GRAPH_DIR=graphs
GRAPH=$1
GRAPH_PATH=$GRAPH_DIR/$GRAPH

hadoop dfs -rmr $GRAPH.sssp.results

time hadoop jar giraph-examples/target/giraph-examples-1.1.0-SNAPSHOT-for-hadoop-1.0.2-jar-with-dependencies.jar org.apache.giraph.GiraphRunner org.apache.giraph.examples.SimpleShortestPathsComputation -vif org.apache.giraph.io.formats.JsonLongDoubleFloatDoubleVertexInputFormat -vip $GRAPH_PATH.json -of org.apache.giraph.io.formats.IdWithValueTextOutputFormat -op $GRAPH.sssp.results -w $WORKERS
