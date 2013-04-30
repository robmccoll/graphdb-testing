#!/bin/bash
CPUMODEL=$(grep -m 1 "model name" /proc/cpuinfo | sed -e 's/.*: //')
SOCKETCOUNT=$(cat /proc/cpuinfo | grep "physical id" | sort | uniq | wc -l)
MEMTOTAL=$(grep -m 1 "MemTotal" /proc/meminfo | sed -e 's/.*:[ \t]*//')
echo "\"sysconfig\" : {
  \"cpu\" : \"$CPUMODEL\",
  \"sockets\": \"$SOCKETCOUNT\",
  \"mem\" : \"$MEMTOTAL\"
}"
