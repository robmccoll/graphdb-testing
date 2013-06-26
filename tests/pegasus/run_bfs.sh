# Program : run_ccmpt.sh
# Description : Run HCC, a connected component algorithm on hadoop.

which hadoop > /dev/null
status=$?
if test $status -ne 0 ; then
	echo ""
	echo "Hadoop is not installed in the system."
	echo "Please install Hadoop and make sure the hadoop binary is accessible."
	exit 127
fi


if [ $# -ne 3 ]; then
	 echo 1>&2 Usage: $0 [#_of_nodes] [#_of_reducers] [HDFS edge_file_path]
	 echo 1>&2 [#_of_nodes] : number of nodes in the graph
	 echo 1>&2 [#_of_reducers] : number of reducers to use in hadoop
	 echo 1>&2 [HDFS edge_file_path] : HDFS directory where edge file is located
	 echo 1>&2    ex: $0 6 3 cc_edge
	 exit 127
fi

rm -rf sssp_output_temp

hadoop dfs -rmr sssp_csr
hadoop dfs -rmr sssp_temp
hadoop dfs -rmr sssp_next
hadoop dfs -rmr sssp_output
hadoop dfs -rmr sssp_summaryout

time hadoop jar pegasus-2.0.jar pegasus.SSSP $3 sssp_csr sssp_temp sssp_next sssp_output $1 $2 nosym

rm -rf sssp_output_temp
