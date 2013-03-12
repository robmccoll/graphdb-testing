package edu.gt.gtri.orientdbtesting;

import java.io.*;
import java.lang.Integer;
import java.util.LinkedList;
import java.util.HashMap;

import com.orientechnologies.orient.core.db.graph.OGraphDatabase;
import com.orientechnologies.orient.core.record.impl.ODocument;
import com.orientechnologies.orient.core.db.record.OIdentifiable;

public class App 
{
  public static void main( String[] args )
  {
    long startTime, endTime;
    long nv = 0, ne = 0;
    long [] off = null, ind = null, wgt = null;
    long na = 0;
    long [] actions = null;
    ODocument [] vertices = null;

    System.out.println( "Up and running..." );

    if(args.length < 2) {
      System.out.printf("Usage: make ARGS=\"graph.file actions.file\"\n");
      System.exit(-1);
    }

    OGraphDatabase graph = new OGraphDatabase("memory:test");
    try {
      graph.create();

      System.out.printf("Reading graph from disk... %s\n", args[0]);
      try {
	FileInputStream graphIn = new FileInputStream(args[0]);

	long endianCheck = readLittleLong(graphIn);
	if(endianCheck != 0x1234ABCD) {
	  System.err.println("Endianness failure " + String.format("%02X", endianCheck));
	  System.exit(-1);
	}

	nv = readLittleLong(graphIn);
	ne = readLittleLong(graphIn);

	System.out.println("\tGraph size nv " + nv + " ne " + ne);

	off = new long[(int)(nv+1)];
	ind = new long[(int)(ne)];
	wgt = new long[(int)(ne)];

	for(int v = 0; v < nv+1; v++) {
	  off[v] = readLittleLong(graphIn);
	}
	for(int e = 0; e < ne; e++) {
	  ind[e] = readLittleLong(graphIn);
	}
	for(int e = 0; e < ne; e++) {
	  wgt[e] = readLittleLong(graphIn);
	}

      } catch (Exception e) {
	System.err.println("Error reading graph.");
	System.exit(-1);
      }
      System.out.println("\tDone\n");

      System.out.println("Loading graph into DB...");
      {
	startTime = System.nanoTime();

	vertices = new ODocument[(int)nv];

	for(long v = 0; v < nv; v++) {
	  vertices[(int)v] = graph.createVertex().field("id", v);
	}
	for(int v = 0; v < nv; v++) {
	  for(int e = (int)off[v]; e < off[v+1]; e++) {
	    graph.createEdge(vertices[v], vertices[(int)ind[e]]).field("weight", 1).save();
	  }
	}

	endTime = System.nanoTime();
      }
      System.out.println("\tDone... " + (((double)(endTime - startTime))/1e9) + "\n");

      System.out.println("Shiloach-Vishkin...");
      {
	startTime = System.nanoTime();

	int [] componentID = new int[(int)nv];
	for(int v = 0; v < nv; v++) {
	  componentID[v] = v;
	}
	
	while(true) {
	  boolean changed = false;
	  for(ODocument edge : graph.browseEdges()) {
	    long srcId = graph.getInVertex(edge).field("id");
	    long dstId = graph.getOutVertex(edge).field("id");

	    if(componentID[(int)srcId] < componentID[(int)dstId]) {
	      componentID[(int)srcId] = componentID[(int)dstId];
	      changed = true;
	    }
	  }

	  if(!changed)
	    break;

	  for(int v = 0; v < nv; v++) {
	    while(componentID[v] != componentID[componentID[v]]) {
	      componentID[v] = componentID[componentID[v]];
	    }
	  }
	}

	endTime = System.nanoTime();

	long count = 0;
	for(int v = 0; v < nv; v++) {
	  if(componentID[v] == v)
	    count++;
	}
	System.out.println("\tComponents... " + count);
      }
      System.out.println("\tDone... " + (((double)(endTime - startTime))/1e9) + "\n");

      System.out.println("BFS...");
      {
	startTime = System.nanoTime();

	HashMap<Long, Boolean> found = new HashMap<Long, Boolean>((int)nv);
	HashMap<Long, Integer> distance = new HashMap<Long, Integer>((int)nv);
	
	for(long v = 1; v < nv; v++) {
	  found.put(v, false);
	  distance.put(v, -1);
	}
	found.put((long)0, true);
	distance.put((long)0, 0);

	LinkedList<ODocument> q = new LinkedList<ODocument>();
	q.add(vertices[0]);
	
	while(!q.isEmpty()) {
	  ODocument v = q.poll();
	  long vid = (Long)v.field("id");
	  for(OIdentifiable edge : graph.getOutEdges(v)) {
	    ODocument endVertex = graph.getInVertex(edge.getRecord());
	    long eid = (Long)endVertex.field("id");
	    if(!found.get(eid)) {
	      found.put(eid, true);
	      distance.put(eid, (distance.get(vid) + 1));
	      q.offer(endVertex);
	    }
	  }
	}

	endTime = System.nanoTime();

	int depth = 0;
	for(int d : distance.values()) {
	  if(d > depth)
	    depth = d;
	}
	System.out.println("\tDepth... " + depth);
      }
      System.out.println("\tDone... " + (((double)(endTime - startTime))/1e9) + "\n");

      System.out.println("Page Rank...");
      {
	double epsilon = 1e-8;
	double dampingfactor = 0.85;
	int maxiter = 100;

	startTime = System.nanoTime();
	HashMap<Long, Double> pr = new HashMap<Long, Double>((int)nv);
	HashMap<Long, Double> prTmp = new HashMap<Long, Double>((int)nv);

	for(long v = 0; v < nv; v++) {
	  pr.put(v, 1.0/((double)nv));
	}

	int iter = maxiter;
	double delta = 1;

	while(delta > epsilon && iter > 0) {
	  for(ODocument v : graph.browseVertices()) {
	    double myPrTmp = 0;
	    for(OIdentifiable edge : graph.getOutEdges(v)) {
	      ODocument u = graph.getInVertex(edge.getRecord());
	      myPrTmp += pr.get((Long)u.field("id")) / ((double)(graph.getOutEdges(u).size()));
	    }
	    prTmp.put((Long)v.field("id"), myPrTmp);
	  }

	  for(ODocument v : graph.browseVertices()) {
	    long id = (Long)v.field("id");
	    prTmp.put(id, prTmp.get(id) * dampingfactor + ((1.0-dampingfactor) / ((double)nv)));
	  }

	  delta = 0;
	  for(ODocument v : graph.browseVertices()) {
	    long id = (Long)v.field("id");
	    double mydelta = prTmp.get(id) - pr.get(id);

	    if(mydelta < 0)
	      mydelta = -mydelta;

	    delta += mydelta;
	    pr.put(id, prTmp.get(id));
	  }

	  //System.out.println("\tIteration " + (maxiter - iter) + " delta " + delta);
	  iter--;
	}

	endTime = System.nanoTime();

	System.out.println("\tIterations... " + (maxiter - iter));
      }
      System.out.println("\tDone... " + (((double)(endTime - startTime))/1e9) + "\n");

      System.out.printf("Reading actions from disk... %s\n", args[1]);
      try {
	FileInputStream actionIn = new FileInputStream(args[1]);

	long endianCheck = readLittleLong(actionIn);
	if(endianCheck != 0x1234ABCD) {
	  System.err.println("Endianness failure " + String.format("%02X", endianCheck));
	  System.exit(-1);
	}

	na = readLittleLong(actionIn);

	System.out.println("\tNum actions " + na);

	actions = new long[(int)na*2];

	for(int a = 0; a < na; a++) {
	  actions[2*a] = readLittleLong(actionIn);
	  actions[2*a+1] = readLittleLong(actionIn);
	}
      } catch (Exception e) {
	System.err.println("Error reading actions.");
	System.exit(-1);
      }
      System.out.println("\tDone\n");

      System.out.println("Updates...");
      {
	startTime = System.nanoTime();

	for(int a = 0; a < na; a++) {
	  long i = actions[2*a];
	  long j = actions[2*a+1];
	  if(i < 0) {
	    i = ~i;
	    j = ~j;
	    for(OIdentifiable edge : graph.getEdgesBetweenVertexes(vertices[(int)i], vertices[(int)j])) {
	      graph.removeEdge(edge.getRecord());
	      break;
	    }
	    for(OIdentifiable edge : graph.getEdgesBetweenVertexes(vertices[(int)j], vertices[(int)i])) {
	      graph.removeEdge(edge.getRecord());
	      break;
	    }
	  } else {
	    boolean handled = false;
	    for(OIdentifiable edge : graph.getEdgesBetweenVertexes(vertices[(int)i], vertices[(int)j])) {
	      ODocument edgedoc = edge.getRecord();
	      edgedoc.field("weight", ((Integer)edgedoc.field("weight")) + 1);
	      break;
	    }
	    if(!handled) {
	      graph.createEdge(vertices[(int)i], vertices[(int)j]).field("weight", 1);
	    }
	    handled = false;
	    for(OIdentifiable edge : graph.getEdgesBetweenVertexes(vertices[(int)j], vertices[(int)i])) {
	      ODocument edgedoc = edge.getRecord();
	      edgedoc.field("weight", ((Integer)edgedoc.field("weight")) + 1);
	      break;
	    }
	    if(!handled) {
	      graph.createEdge(vertices[(int)j], vertices[(int)i]).field("weight", 1);
	    }
	  }
	}

	endTime = System.nanoTime();
      }
      System.out.println("\tDone... " + (((double)(endTime - startTime))/1e9) + "\n");
    } finally {
      graph.close();
    }
  }

  public static long readLittleLong(FileInputStream in) throws IOException {
    byte[] b = new byte[8];
    try {
      in.read(b);
    } catch (Exception e) {
      System.err.println("Error reading little long.\n");
    }
    return  (((long)b[7] & 0xffL) << 56) +
	    (((long)b[6] & 0xffL) << 48) +
	    (((long)b[5] & 0xffL) << 40) +
	    (((long)b[4] & 0xffL) << 32) +
	    (((long)b[3] & 0xffL) << 24) +
	    (((long)b[2] & 0xffL) << 16) +
	    (((long)b[1] & 0xffL) << 8) +
	     ((long)b[0] & 0xffL);
  }
}
