#!/bin/bash

WORKERS=$3

GRAPH_DIR=graphs
GRAPH=$1
GRAPH_PATH=$GRAPH_DIR/$GRAPH

hadoop dfs -rmr $GRAPH.insert.results

#time hadoop jar giraph-examples/target/giraph-examples-1.1.0-SNAPSHOT-for-hadoop-1.0.2-jar-with-dependencies.jar org.apache.giraph.GiraphRunner org.apache.giraph.examples.InsertRemoveEdges -vif org.apache.giraph.io.formats.IntIntNullTextInputFormat -vip $GRAPH_PATH.intint -of org.apache.giraph.io.formats.IdWithValueTextOutputFormat -op $GRAPH.insert.results -w $WORKERS -ca InsertRemoveEdges.actionFile=$2

time hadoop jar giraph-examples/target/giraph-examples-1.1.0-SNAPSHOT-for-hadoop-1.0.2-jar-with-dependencies.jar org.apache.giraph.GiraphRunner org.apache.giraph.examples.InsertRemoveEdges\$InsertRemoveEdgesComputation -vif org.apache.giraph.io.formats.IntIntNullTextInputFormat -vip $GRAPH_PATH.intint -of org.apache.giraph.io.formats.IdWithValueTextOutputFormat -op $GRAPH.insert.results -w $WORKERS -ca InsertRemoveEdges.actionFile=$2
