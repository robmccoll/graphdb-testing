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
File: SSSP.java
- HCC: Do breadth-first search from vertex 0
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


class SSSPResultInfo
{
  public int white;
  public int gray;
  public int black;
};

public class SSSP extends Configured implements Tool 
{
  public static int MAX_ITERATIONS = 2048;
  public static int white_nodes[] = new int[MAX_ITERATIONS];
  public static int gray_nodes[] = new int[MAX_ITERATIONS];
  public static int black_nodes[] = new int[MAX_ITERATIONS];
  static int iter_counter = 0;

  //////////////////////////////////////////////////////////////////////
  // STAGE 1: Form adjacency representation of the graph.
  //	- Format is "Vtx Clr Dist Neighbors"
  //	- Clr: Unvisited = 0, Processing = 1, Done = 2
  //////////////////////////////////////////////////////////////////////
  public static class MapStage1 extends MapReduceBase	implements Mapper<LongWritable, Text, IntWritable, Text>
  {
    //private final IntWritable from_node_int = new IntWritable();
    //private final IntWritable to_node_int = new IntWritable();
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

      //from_node_int.set(Integer.parseInt(line[0]));
      //to_node_int.set(Integer.parseInt(line[1]));
      from_node_int = Integer.parseInt(line[0]);
      to_node_int = Integer.parseInt(line[1]);
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

      if (key.get() == 0) {    // source vertex for the search
	vertex_id_str = "" + src + "\t" + 1 + "\t" + 0;
      } else {
	vertex_id_str = "" + src + "\t" + 0 + "\t" + Integer.MAX_VALUE;
      }
      
      while (values.hasNext()) {
	Text from_cur_node = values.next();
	String line = from_cur_node.toString();
	int neighbor = Integer.parseInt(line);
	vertex_id_str = vertex_id_str + "\t" + neighbor;
      }

      output.collect(key, new Text(vertex_id_str));
    }
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // STAGE 2: Process color = 1, producing next frontier
  ////////////////////////////////////////////////////////////////////////////////////////////////
  public static class MapStage2 extends MapReduceBase implements Mapper<LongWritable, Text, IntWritable, Text>
  {
    // Identity mapper
    public void map (final LongWritable key, final Text value, final OutputCollector<IntWritable, Text> output, final Reporter reporter) throws IOException
    {
      final String[] line = value.toString().split("\t");
      String vertex_id_str = "";

      int src = Integer.parseInt(line[1]);
      int color = Integer.parseInt(line[2]);
      int distance = Integer.parseInt(line[3]);
      //System.out.println("   " + key.get() + "\t" + src + "\t" + color + "\t" + distance);

      int degree = line.length-4;

      if (color == 0) {
	vertex_id_str = "" + src + "\t" + 0 + "\t" + distance;
	for (int i = 0; i < degree; i++) {
	  vertex_id_str = vertex_id_str + "\t" + Integer.parseInt(line[4+i]);
	}
	output.collect(new IntWritable(src), new Text(vertex_id_str));
      }
      else if (color == 1) {
	for (int i = 0; i < degree; i++) {
	  int neighbor = Integer.parseInt(line[4+i]);
	  if (neighbor != src) {
	    vertex_id_str = "" + neighbor + "\t" + 1 + "\t" + (distance + 1);
	    output.collect(new IntWritable(neighbor), new Text(vertex_id_str));
	  }
	}
	vertex_id_str = "" + src + "\t" + 2 + "\t" + distance;
	output.collect(new IntWritable(src), new Text(vertex_id_str));
      } else if (color == 2) {
	vertex_id_str = "" + src + "\t" + 2 + "\t" + distance;
	output.collect(new IntWritable(src), new Text(vertex_id_str));
      }
    }
  }

  public static class RedStage2 extends MapReduceBase	implements Reducer<IntWritable, Text, IntWritable, Text>
  {
    public void reduce (final IntWritable key, final Iterator<Text> values, final OutputCollector<IntWritable, Text> output, final Reporter reporter) throws IOException
    {
      int src = key.get();
      String vertex_id_str = "";
      String adjacency_list = "";

      int newColor = 0;
      int newDistance = Integer.MAX_VALUE;
      
      while (values.hasNext()) {
	Text cur_value_text = values.next();
	final String[] line = cur_value_text.toString().split("\t");

	int color = Integer.parseInt(line[1]);
	if (color > newColor)
	  newColor = color;

	int distance = Integer.parseInt(line[2]);
	if (distance < newDistance)
	  newDistance = distance;

	int degree = line.length-3;
	for (int i = 0; i < degree; i++) {
	  adjacency_list = adjacency_list + "\t" + Integer.parseInt(line[3+i]);
	}
      }

      vertex_id_str = "" + src + "\t" + newColor + "\t" + newDistance + adjacency_list;
      output.collect(key, new Text(vertex_id_str));
      
/*
      while (values.hasNext())
      {
	output.collect(key, values.next());
      }
      */
    }
  }


  //////////////////////////////////////////////////////////////////////
  // STAGE 3: Calculate number of vertices on the next frontier.
  //////////////////////////////////////////////////////////////////////
  public static class	MapStage3 extends MapReduceBase implements Mapper<LongWritable, Text, IntWritable, IntWritable>
  {
    // output : f n		( n : # of node whose component didn't change)
    //          i m		( m : # of node whose component changed)
    public void map (final LongWritable key, final Text value, final OutputCollector<IntWritable, IntWritable> output, final Reporter reporter) throws IOException
    {
      if (value.toString().startsWith("#"))
	return;

      final String[] line = value.toString().split("\t");
      int color = Integer.parseInt(line[2]);

      output.collect(new IntWritable(color), new IntWritable(1) );
    }
  }

  public static class	RedStage3 extends MapReduceBase	implements Reducer<IntWritable, IntWritable, Text, Text>
  {
    public void reduce (final IntWritable key, final Iterator<IntWritable> values, final OutputCollector<Text, Text> output, final Reporter reporter) throws IOException
    {
      int sum = 0;
      int color = key.get();

      while (values.hasNext()) {
	final IntWritable line = values.next();
	sum += 1;
      }

      output.collect(new Text(Integer.toString(color)), new Text(Integer.toString(sum)) );
    }
  }


  //////////////////////////////////////////////////////////////////////
  // command line interface
  //////////////////////////////////////////////////////////////////////
  protected Path edge_path = null;
  protected Path csr_path = null;
  protected Path tempbm_path = null;
  protected Path nextbm_path = null;
  protected Path output_path = null;
  protected Path summaryout_path = null;
  protected String local_output_path;
  protected int number_nodes = 0;
  protected int nreducers = 1;
  protected int cur_iter = 1;
  protected int make_symmetric = 0;		// convert directed graph to undirected graph

  // Main entry point.
  public static void main (final String[] args) throws Exception
  {
    final int result = ToolRunner.run(new Configuration(), new SSSP(), args);

    System.exit(result);
  }

  // Print the command-line usage text.
  protected static int printUsage ()
  {
    System.out.println("sssp <edge_path> <csr_path> <tempbm_path> <nextbm_path> <output_path> <# of nodes> <# of tasks> <makesym or nosym>");

    ToolRunner.printGenericCommandUsage(System.out);

    return -1;
  }

  // submit the map/reduce job.
  public int run (final String[] args) throws Exception
  {
    if( args.length != 8 ) {
      return printUsage();
    }

    edge_path = new Path(args[0]);
    csr_path = new Path(args[1]);
    tempbm_path = new Path(args[2]);
    nextbm_path = new Path(args[3]);
    output_path = new Path(args[4]);
    summaryout_path = new Path("concmpt_summaryout");
    number_nodes = Integer.parseInt(args[5]);
    nreducers = Integer.parseInt(args[6]);

    if( args[7].compareTo("makesym") == 0 )
      make_symmetric = 1;
    else
      make_symmetric = 0;

    System.out.println("\n-----===[PEGASUS: A Peta-Scale Graph Mining System]===-----\n");
    System.out.println("[PEGASUS] Computing breadth-first search. Edge path = " + args[0] + ", Reducers = " + nreducers );

    local_output_path = args[4] + "_temp";

    JobClient.runJob(configStage1());
    FileSystem.get(getConf()).rename(csr_path, tempbm_path);

    // Iteratively calculate neighborhood function. 
    for (int i = cur_iter; i < MAX_ITERATIONS; i++) {
      cur_iter++;

      JobClient.runJob(configStage2());
      JobClient.runJob(configStage3());

      FileUtil.fullyDelete( FileSystem.getLocal(getConf()), new Path(local_output_path));

      final FileSystem fs = FileSystem.get(getConf());

      // copy neighborhood information from HDFS to local disk, and read it!
      String new_path = local_output_path + "/" + i;
      fs.copyToLocalFile(output_path, new Path(new_path) ) ;
      SSSPResultInfo ri = readIterationOutput(new_path);

      white_nodes[iter_counter] = ri.white;
      gray_nodes[iter_counter] = ri.gray;
      black_nodes[iter_counter] = ri.black;

      iter_counter++;

      System.out.println("Hop " + i + " : white = " + ri.white + ", gray = " + ri.gray + ", black = " + ri.black);

      // Stop when the minimum neighborhood doesn't change
      if( ri.gray == 0 ) {
	System.out.println("All vertices have been reached. Finishing...");
	//fs.delete(csr_path);
	fs.delete(tempbm_path);
	fs.delete(output_path);
	fs.rename(nextbm_path, output_path);

	break;
      }

      // rotate directory
      //fs.delete(csr_path);
      fs.delete(tempbm_path);
      fs.delete(output_path);
      fs.rename(nextbm_path, tempbm_path);

    }

 //   FileUtil.fullyDelete( FileSystem.getLocal(getConf()), new Path(local_output_path));

    // calculate summary information using an additional stage
    //System.out.println("Summarizing connected components information...");
    //JobClient.runJob(configStage4());

    // finishing.
    System.out.println("\n[PEGASUS] Breadth-first search computed.");
    System.out.println("[PEGASUS] Total Iteration = " + iter_counter);
    System.out.println("[PEGASUS] BFS distance labels are saved in the HDFS sssp_output as\n\"vertex	vertex	number	distance\" format.\n");

    return 0; 
  }

  // read neighborhood number after each iteration.
  public static SSSPResultInfo readIterationOutput(String new_path) throws Exception
  {
    SSSPResultInfo ri = new SSSPResultInfo();
    ri.white = ri.gray = ri.black = 0;
    String output_path = new_path + "/part-00000";
    String file_line = "";
    //System.out.println("Reading " + output_path + "...");

    try {
      BufferedReader in = new BufferedReader(	new InputStreamReader(new FileInputStream( output_path ), "UTF8"));

      // Read first line
      file_line = in.readLine();

      // Read through file one line at time. Print line # and line
      while (file_line != null){
	final String[] line = file_line.split("\t");

	//System.out.println("\t" + Integer.parseInt(line[0]) + ":\t" + Integer.parseInt(line[1]));
	int color = Integer.parseInt(line[0]);
	int number = Integer.parseInt(line[1]);
	if (color == 0)
	  ri.white = number;
	else if (color == 1)
	  ri.gray = number;
	else if (color == 2)
	  ri.black = number;
	else
	  System.out.println("ERROR counting colors");

	file_line = in.readLine();
      }

      in.close();
    } catch (IOException e) {
      e.printStackTrace();
    }

    return ri;//result;
  }

  // Configure stage1
  protected JobConf configStage1() throws Exception
  {
    final JobConf conf = new JobConf(getConf(), SSSP.class);
    conf.set("number_nodes", "" + number_nodes);
    conf.set("cur_iter", "" + cur_iter);
    conf.set("make_symmetric", "" + make_symmetric);
    conf.setJobName("SSSP_Stage1");

    conf.setMapperClass(MapStage1.class);
    conf.setReducerClass(RedStage1.class);

    FileInputFormat.setInputPaths(conf, edge_path);  
    FileOutputFormat.setOutputPath(conf, csr_path);  

    conf.setNumReduceTasks( nreducers );

    conf.setOutputKeyClass(IntWritable.class);
    conf.setOutputValueClass(Text.class);

    return conf;
  }

  // Configure stage2
  protected JobConf configStage2 () throws Exception
  {
    final JobConf conf = new JobConf(getConf(), SSSP.class);
    conf.set("number_nodes", "" + number_nodes);
    conf.set("cur_iter", "" + cur_iter);
    conf.set("make_symmetric", "" + make_symmetric);
    conf.setJobName("SSSP_Stage2");

    conf.setMapperClass(MapStage2.class);        
    conf.setReducerClass(RedStage2.class);
    //conf.setCombinerClass(CombinerStage2.class);

    FileInputFormat.setInputPaths(conf, tempbm_path);  
    FileOutputFormat.setOutputPath(conf, nextbm_path);  

    conf.setNumReduceTasks( nreducers );

    conf.setOutputKeyClass(IntWritable.class);
    conf.setOutputValueClass(Text.class);

    return conf;
  }

  // Configure stage3
  protected JobConf configStage3 () throws Exception
  {
    final JobConf conf = new JobConf(getConf(), ConCmpt.class);
    conf.set("number_nodes", "" + number_nodes);
    conf.setJobName("SSSP_Stage3");

    conf.setMapperClass(MapStage3.class);        
    conf.setReducerClass(RedStage3.class);
    //conf.setCombinerClass(RedStage3.class);

    FileInputFormat.setInputPaths(conf, nextbm_path);  
    FileOutputFormat.setOutputPath(conf, output_path);  

    conf.setNumReduceTasks( 1 );	// This is necessary.

    conf.setOutputKeyClass(IntWritable.class);
    conf.setOutputValueClass(IntWritable.class);

    return conf;
  }

}

