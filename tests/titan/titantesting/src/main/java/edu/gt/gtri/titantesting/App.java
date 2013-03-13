package edu.gt.gtri.titantesting;

import java.io.*;
import java.lang.Integer;
import java.util.LinkedList;
import java.util.HashMap;

import com.thinkaurelius.titan.core.TitanGraph;
import com.thinkaurelius.titan.core.TitanFactory;
import com.tinkerpop.blueprints.Vertex;
import com.tinkerpop.blueprints.Edge;
import com.tinkerpop.blueprints.Direction;
import com.tinkerpop.blueprints.Graph;

public class App 
{
  public static void main( String[] args )
  {
    long startTime, endTime;
    long nv = 0, ne = 0;
    long [] off = null, ind = null, wgt = null;
    long na = 0;
    long [] actions = null;
    TitanGraph g = null;
    Vertex [] vertices = null;

    System.out.println( "Up and running..." );

    if(args.length < 2) {
      System.out.printf("Usage: make ARGS=\"graph.file actions.file\"\n");
      System.exit(-1);
    }


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

    vertices = new Vertex[(int)nv];

    System.out.println("Loading graph into DB...");
    {
      startTime = System.nanoTime();
      g = TitanFactory.openInMemoryGraph();

      for(int v = 0; v < nv; v++) {
	vertices[(int)v] = g.addVertex(null);
	vertices[(int)v].setProperty("vid", v);
      }
      for(int v = 0; v < nv; v++) {
	Vertex src = vertices[(int)v];
	for(int e = (int)off[v]; e < off[v+1]; e++) {
	  g.addEdge(null, src, vertices[(int)ind[e]], "edge").setProperty("weight", 1);
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
	for(Edge e : g.getEdges()) {
	  int src = (Integer)e.getVertex(Direction.OUT).getProperty("vid");
	  int dst = (Integer)e.getVertex(Direction.IN).getProperty("vid");

	  if(componentID[src] < componentID[dst]) {
	    componentID[src] = componentID[dst];
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

      HashMap<Vertex, Boolean> found = new HashMap<Vertex, Boolean>((int)nv);
      HashMap<Vertex, Integer> distance = new HashMap<Vertex, Integer>((int)nv);
      LinkedList<Vertex> q = new LinkedList<Vertex>();

      for(Vertex v : g.getVertices()) {
	found.put(v, false);
	distance.put(v, -1);
      }

      found.put(vertices[0], true);
      distance.put(vertices[0], 0);
      q.add(vertices[0]);

      while(!q.isEmpty()) {
	Vertex v = q.poll();
	for(Vertex u : v.query().labels("edge").vertices()) {
	  if(!found.get(u)) {
	    found.put(u, true);
	    distance.put(u, distance.get(v) + 1);
	    q.offer(u);
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
      HashMap<Vertex, Double> pr = new HashMap<Vertex, Double>((int)nv);
      HashMap<Vertex, Double> prTmp = new HashMap<Vertex, Double>((int)nv);

      for(Vertex v : g.getVertices()) {
	pr.put(v, 1.0/((double)nv));
      }

      int iter = maxiter;
      double delta = 1;

      while(delta > epsilon && iter > 0) {
	for(Vertex v : g.getVertices()) {
	  double myPrTmp = 0;
	  for(Vertex u : v.query().labels("edge").vertices()) {
	    myPrTmp += pr.get(u) / ((double)u.query().labels("edge").count());
	  }
	  prTmp.put(v, myPrTmp);
	}

	for(Vertex v : g.getVertices()) {
	  prTmp.put(v, prTmp.get(v) * dampingfactor + ((1.0-dampingfactor) / ((double)nv)));
	}

	delta = 0;
	for(Vertex v : g.getVertices()) {
	  double mydelta = prTmp.get(v) - pr.get(v);

	  if(mydelta < 0)
	    mydelta = -mydelta;

	  delta += mydelta;
	  pr.put(v, prTmp.get(v));
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
	  i =	~i;
	  j =	~j;
	  for(Edge e : vertices[(int)i].query().labels("edge").edges()) {
	    if(j == (Integer)e.getVertex(Direction.OUT).getProperty("vid")) {
	      g.removeEdge(e);
	      break;
	    }
	  }
	  for(Edge e : vertices[(int)j].query().labels("edge").edges()) {
	    if(i == (Integer)e.getVertex(Direction.OUT).getProperty("vid")) {
	      g.removeEdge(e);
	      break;
	    }
	  }
	} else {
	  boolean handled = false;
	  for(Edge e : vertices[(int)i].query().labels("edge").edges()) {
	    if(j == (Integer)e.getVertex(Direction.OUT).getProperty("vid")) {
	      e.setProperty("weight", 1 + (Integer)e.getProperty("weight"));
	      handled = true;
	      break;
	    }
	  }
	  if(!handled) {
	    g.addEdge(null, vertices[(int)i], vertices[(int)j], "edge").setProperty("weight", 1);
	  }
	  handled = false;
	  for(Edge e : vertices[(int)j].query().labels("edge").edges()) {
	    if(i == (Integer)e.getVertex(Direction.OUT).getProperty("vid")) {
	      e.setProperty("weight", 1 + (Integer)e.getProperty("weight"));
	      handled = true;
	      break;
	    }
	  }
	  if(!handled) {
	    g.addEdge(null, vertices[(int)j], vertices[(int)i], "edge").setProperty("weight", 1);
	  }
	}
      }

      endTime = System.nanoTime();
    }
    System.out.println("\tDone... " + (((double)(endTime - startTime))/1e9) + "\n");
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
