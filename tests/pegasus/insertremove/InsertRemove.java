/***********************************************************************
PEGASUS: Peta-Scale Graph Mining System
Authors: U Kang, Duen Horng Chau, and Christos Faloutsos

This software is licensed under Apache License, Version 2.0 (the  "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
-------------------------------------------------------------------------
File: InsertRemove.java
Version: 2.0
 ***********************************************************************/

package pegasus;

import java.io.*;
import java.util.*;
import java.text.*;
import org.apache.hadoop.conf.*;
import org.apache.hadoop.fs.*;
import org.apache.hadoop.io.*;
import org.apache.hadoop.mapred.*;
import org.apache.hadoop.util.*;


public class InsertRemove extends Configured implements Tool 
{

  //////////////////////////////////////////////////////////////////////
  // STAGE 1: Form adjacency representation of the graph.
  //////////////////////////////////////////////////////////////////////
  public static class MapStage1 extends MapReduceBase	implements Mapper<LongWritable, Text, IntWritable, Text>
  {
    int from_node_int = 0;
    int to_node_int = 0; 
    int make_symmetric = 0;

    public void configure(JobConf job) {
      make_symmetric = Integer.parseInt(job.get("make_symmetric"));

      System.out.println("MapStage1 : make_symmetric = " + make_symmetric);
    }

    public void map (final LongWritable key, final Text value, final OutputCollector<IntWritable, Text> output, final Reporter reporter) throws IOException
    {
      String line_text = value.toString();
      if (line_text.startsWith("#"))				// ignore comments in the edge file
	return;

      final String[] line = line_text.split("\t");
      if( line.length < 2 )
	return;

      from_node_int = Integer.parseInt(line[0]);
      to_node_int = Integer.parseInt(line[1]);
      if (from_node_int >= 0)
	output.collect(new IntWritable(from_node_int), new Text(Integer.toString(to_node_int)));

      if (make_symmetric == 1) {
	output.collect(new IntWritable(to_node_int), new Text(Integer.toString(from_node_int)));
      }
    }
  }

  public static class	RedStage1 extends MapReduceBase implements Reducer<IntWritable, Text, IntWritable, Text>
  {
    int number_nodes = 0;

    public void configure(JobConf job) {
      number_nodes = Integer.parseInt(job.get("number_nodes"));

      System.out.println("RedStage1 : configure is called. number_nodes = " + number_nodes );
    }

    public void reduce (final IntWritable key, final Iterator<Text> values, OutputCollector<IntWritable, Text> output, final Reporter reporter) throws IOException
    {
      String vertex_id_str = "";
      int src = key.get();

      vertex_id_str = "" + src;
      
      while (values.hasNext()) {
	Text from_cur_node = values.next();
	String line = from_cur_node.toString();
	int neighbor = Integer.parseInt(line);
	vertex_id_str = vertex_id_str + "\t" + neighbor;
      }

      output.collect(key, new Text(vertex_id_str));
    }
  }




  //////////////////////////////////////////////////////////////////////
  // command line interface
  //////////////////////////////////////////////////////////////////////
  protected Path edge_path = null;
  protected Path action_path = null;
  protected Path csr_path = null;
  protected int number_nodes = 0;
  protected int nreducers = 1;
  protected int make_symmetric = 0;		// convert directed graph to undirected graph

  // Main entry point.
  public static void main (final String[] args) throws Exception
  {
    final int result = ToolRunner.run(new Configuration(), new InsertRemove(), args);

    System.exit(result);
  }

  // Print the command-line usage text.
  protected static int printUsage ()
  {
    System.out.println("insertremove <edge_path> <action_path> <csr_path> <# of nodes> <# of tasks> <makesym or nosym>");

    ToolRunner.printGenericCommandUsage(System.out);

    return -1;
  }

  // submit the map/reduce job.
  public int run (final String[] args) throws Exception
  {
    if( args.length != 6 ) {
      return printUsage();
    }

    edge_path = new Path(args[0]);
    action_path = new Path(args[1]);
    csr_path = new Path(args[2]);
    number_nodes = Integer.parseInt(args[3]);
    nreducers = Integer.parseInt(args[4]);

    if( args[5].compareTo("makesym") == 0 )
      make_symmetric = 1;
    else
      make_symmetric = 0;

    System.out.println("\n-----===[PEGASUS: A Peta-Scale Graph Mining System]===-----\n");
    System.out.println("[PEGASUS] Running edge insertions and building the graph. Edge path = " + args[0] + ", Action path = " + args[1] + ", Reducers = " + nreducers );

    JobClient.runJob(configStage1());


    System.out.println("\n[PEGASUS] All insertions have been processed.");
    System.out.println("[PEGASUS] New compressed sparse row graph is stored in HDFS csr_graph.\n");

    return 0; 
  }


  // Configure stage1
  protected JobConf configStage1() throws Exception
  {
    final JobConf conf = new JobConf(getConf(), InsertRemove.class);
    conf.set("number_nodes", "" + number_nodes);
    conf.set("make_symmetric", "" + make_symmetric);
    conf.setJobName("InsertRemove_Stage1");

    conf.setMapperClass(MapStage1.class);
    conf.setReducerClass(RedStage1.class);

    FileInputFormat.setInputPaths(conf, edge_path, action_path);  
    FileOutputFormat.setOutputPath(conf, csr_path);  

    conf.setNumReduceTasks( nreducers );

    conf.setOutputKeyClass(IntWritable.class);
    conf.setOutputValueClass(Text.class);

    return conf;
  }


}

