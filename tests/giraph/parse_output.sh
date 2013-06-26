#!/bin/bash

GRAPH_FILE=$1
NAME=giraph
VER=$2
NV=0
NE=0
NA=0

if [[ $GRAPH_FILE == "tiny.g" ]]
then
  NV=1024
  NE=14408
  NA=100000
elif [[ $GRAPH_FILE == "small.g" ]]
then
  NV=32768
  NE=504364
  NA=100000
elif [[ $GRAPH_FILE == "medium.g" ]]
then
  NV=1048576
  NE=16571375
  NA=1000000
elif [[ $GRAPH_FILE == "large.g" ]]
then
  NV=16777216
  NE=267069269
  NA=1000000
fi



echo "{"
echo "\"type\":\"$NAME\","
echo "\"nv\":$NV,"
echo "\"ne\":$NE,"

MILLISECONDS=$(cat "$GRAPH_FILE.cc.$VER.err" | grep milliseconds | grep "Input superstep" | cut -d= -f2 | awk '{s+=$1} END {print s}')
MYVAR=`echo "scale=3; $MILLISECONDS/1000" | bc | awk '{printf "%.3f", $0}'`

echo "\"results\": {"
echo "\"build\": {"
echo "\"name\":\"$NAME\","
echo "\"time\":$MYVAR"
echo "},"

MILLISECONDS=$(cat "$GRAPH_FILE.cc.$VER.err" | grep milliseconds | grep Superstep | cut -d= -f2 | awk '{s+=$1} END {print s}')
MYVAR=`echo "scale=3; $MILLISECONDS/1000" | bc | awk '{printf "%.3f", $0}'`

echo "\"sv\": {"
echo "\"name\":\"$NAME\","
echo "\"time\":$MYVAR"
echo "},"

MILLISECONDS=$(cat "$GRAPH_FILE.sssp.$VER.err" | grep milliseconds | grep Superstep | cut -d= -f2 | awk '{s+=$1} END {print s}')
MYVAR=`echo "scale=3; $MILLISECONDS/1000" | bc | awk '{printf "%.3f", $0}'`

echo "\"sssp\": {"
echo "\"name\":\"$NAME\","
echo "\"time\":$MYVAR"
echo "},"

MILLISECONDS=$(cat "$GRAPH_FILE.pr.$VER.err" | grep milliseconds | grep Superstep | cut -d= -f2 | awk '{s+=$1} END {print s}')
MYVAR=`echo "scale=3; $MILLISECONDS/1000" | bc | awk '{printf "%.3f", $0}'`

echo "\"pr\": {"
echo "\"name\":\"$NAME\","
echo "\"time\":$MYVAR"
echo "},"

MILLISECONDS=$(cat "$GRAPH_FILE.insert.$VER.err" | grep milliseconds | grep "Superstep [1-2]" | cut -d= -f2 | awk '{s+=$1} END {print s}')
MYVAR=`echo "scale=3; $NA/($MILLISECONDS/1000)" | bc | awk '{printf "%.3f", $0}'`

echo "\"update\": {"
echo "\"name\":\"$NAME\","
echo "\"time\":$MYVAR"
echo "}"

echo "},"
echo "\"na\":$NA"
#echo "\"mem\":"

echo "}"


#MILLISECONDS=`cat $OUTPUT_FILE | grep milliseconds | grep Superstep | cut -d= -f2 | awk '{s+=$1} END {print s}'`
#SECONDS=`echo "scale=3;$MILLISECONDS/1000" | bc`

