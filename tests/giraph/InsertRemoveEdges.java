/*
 * Licensed to the Apache Software Foundation (ASF) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The ASF licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.apache.giraph.examples;

import java.io.*;
import java.util.*;
import java.text.*;

import org.apache.giraph.graph.BasicComputation;
import org.apache.giraph.conf.StrConfOption;
import org.apache.giraph.edge.Edge;
import org.apache.giraph.edge.EdgeFactory;
import org.apache.giraph.graph.Vertex;
import org.apache.giraph.worker.WorkerContext;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.io.*;
import org.apache.log4j.Logger;

/**
 * Demonstrates the basic capability of inserting and removing edges.
 */
@Algorithm(
    name = "Insert/Remove Edges Benchmark",
    description = "Inserts and removes batches of edges."
    )
public class InsertRemoveEdges //extends BasicComputation<
//IntWritable, IntWritable, NullWritable, IntWritable> {
{
  public static class EdgeActionMessage implements Writable {
    private int source;
    private int dest;
    private int action;

    public EdgeActionMessage() { }

    public EdgeActionMessage(
	int source, int dest, int action) {
      this.source= source;
      this.dest = dest;
      this.action = action;
    }

    @Override
      public void readFields(DataInput input) throws IOException {
	source = input.readInt();
	dest = input.readInt();
	action = input.readInt();
      }

    @Override
      public void write(DataOutput output) throws IOException {
	output.writeInt(source);
	output.writeInt(dest);
	output.writeInt(action);
      }

    @Override
      public String toString() {
	return "(source=" + source + ",dest=" +
	  dest + ",action=" + action + ")";
      }
  }

  public static class InsertRemoveEdgesComputation extends
    BasicComputation<IntWritable, IntWritable, NullWritable,
    EdgeActionMessage> {

      /** The shortest paths id */
      public static final StrConfOption ACTION_FILE =
	new StrConfOption("InsertRemoveEdges.actionFile", "__null__");

      @Override
	public void compute(
	    Vertex<IntWritable, IntWritable, NullWritable> vertex,
	    Iterable<EdgeActionMessage> messages) throws IOException
	{
	  //vertex.setValue(new IntWritable(Integer.MAX_VALUE));
	  //vertex.setValue(new IntWritable(0));

	  /* SUPERSTEP 0 -- read in the action stream */

	  if (getSuperstep() == 0) {
	    vertex.setValue(new IntWritable(0));

	    if (vertex.getId().get() == 1) {
	      int vertexCount = (int) getTotalNumVertices();
	      int edgeCount = (int) getTotalNumEdges();
	      vertex.setValue(new IntWritable(edgeCount));

	      /* Read the list of edges in the action stream and send messages */
	      BufferedReader in = new BufferedReader (new InputStreamReader(new FileInputStream(ACTION_FILE.get(getConf()))));

	      String file_line = in.readLine();

	      Random randGen = new Random();

	      int count = 0;
	      while (file_line != null) {
		final String[] line = file_line.split("\t");

		int randomInt = randGen.nextInt(vertexCount);
		int source = Integer.parseInt(line[0]);
		int dest = Integer.parseInt(line[1]);
		if (source < 0) {
		  sendMessage(new IntWritable(randomInt), new EdgeActionMessage(-source-1, -dest-1, 2));
		}
		else {
		  sendMessage(new IntWritable(randomInt), new EdgeActionMessage(source, dest, 1));
		}
		count++;

		file_line = in.readLine();
	      }

	      in.close();

	    }
	  }

	  /* SUPERSTEP 2 -- assign actions to vertices through messages */

	  if (getSuperstep() == 1) {
	    for (EdgeActionMessage message : messages) {
	      int source = message.source;
	      int dest = message.dest;
	      int action = message.action;
	      sendMessage(new IntWritable(source), new EdgeActionMessage(source, dest, action));
	    }
	  }

	  /* SUPERSTEP 2 -- add new edges */

	  if (getSuperstep() == 2) {
	    if (vertex.getId().get() == 3) {
	      int edgeCount = (int) getTotalNumEdges();
	      vertex.setValue(new IntWritable(edgeCount));
	    }

	    for (EdgeActionMessage message : messages) {
	      int dest = message.dest;
	      int action = message.action;
	      if (action == 2) {
		vertex.removeEdges(new IntWritable(dest));
	      }
	      else {
		Edge<IntWritable, NullWritable> edge = EdgeFactory.create(new IntWritable(dest));
		vertex.addEdge(edge);
	      }
	    }
	    sendMessage(vertex.getId(), new EdgeActionMessage(0, 0, 0));
	  }

	  /* SUPERSTEP 3 -- check edges were added */

	  if (getSuperstep() == 3) {
	    if (vertex.getId().get() == 4) {
	      int edgeCount = (int) getTotalNumEdges();
	      vertex.setValue(new IntWritable(edgeCount));
	    }
	  }

	  vertex.voteToHalt();
	}
    }

}
