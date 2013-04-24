package edu.gt.gtri.neo4jtesting;

import java.io.*;
import java.lang.Integer;
import java.util.LinkedList;
import java.util.HashMap;
import java.lang.management.*;

import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.tooling.GlobalGraphOperations;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Relationship;
import org.neo4j.graphdb.Transaction;
import org.neo4j.graphdb.RelationshipType;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;
import org.neo4j.kernel.EmbeddedGraphDatabase;

/**
 * Hello world!
 *
 */
public class App 
{
  public static GraphDatabaseService graphDb;

  public static enum RelTypes implements RelationshipType {
    EDGE
  }

  public static void R(String result) {
    System.out.println("RSLT: " + result);
  }

  public static void main( String[] args ) {
    long startTime, endTime;
    long nv = 0, ne = 0;
    long [] off = null, ind = null, wgt = null;
    long na = 0;
    long [] actions = null;
    Node [] nodes = null;

    System.out.println( "Up and running..." );

    if(args.length < 2) {
      System.out.printf("Usage: make ARGS=\"graph.file actions.file\"\n");
      System.exit(-1);
    }

    graphDb = new GraphDatabaseFactory().newEmbeddedDatabase("./neo4jdb.db");
    registerShutdownHook(graphDb);


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

      R("{");
      R("\"type\":\"neo4j\",");
      R("\"nv\":" + nv + ",");
      R("\"ne\":" + ne + ",");
      R("\"results\": {");

    System.out.println("Loading graph into DB...");
    {
      startTime = System.nanoTime();

      nodes = new Node[(int)(nv)];

      Transaction tx = graphDb.beginTx();
      try {
	for(int v = 0; v < nv; v++) {
	  nodes[v] = graphDb.createNode();
	  nodes[v].setProperty("id", v);
	}
	for(int v = 0; v < nv; v++) {
	  for(int e = (int)off[v]; e < off[v+1]; e++) {
	    Relationship rel = nodes[v].createRelationshipTo(nodes[(int)ind[e]], RelTypes.EDGE);
	    rel.setProperty("weight", wgt[e]);
	  }
	}

	tx.success();
      } catch (Exception e) {
	tx.failure();
      } finally {
	tx.finish();
      }

      endTime = System.nanoTime();
    }
    double build_time = (((double)(endTime - startTime))/1e9);
    R("\"build\": {");
    R("\"name\":\"neo4j-std\",");
    R("\"time\":" +  build_time);
    R("},");
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
	for(Relationship rel : GlobalGraphOperations.at(graphDb).getAllRelationships()) {
	  Object srcIdObj = rel.getStartNode().getProperty("id", null);
	  Object dstIdObj = rel.getEndNode().getProperty("id", null);

	  if(srcIdObj != null && dstIdObj != null) {
	    int srcId = (Integer)(srcIdObj);
	    int dstId = (Integer)(dstIdObj);
	    if(componentID[srcId] < componentID[dstId]) {
	      componentID[srcId] = componentID[dstId];
	      changed = true;
	    }
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
    double sv_time = (((double)(endTime - startTime))/1e9);
    R("\"sv\": {");
    R("\"name\":\"neo4j-std\",");
    R("\"time\":" +  sv_time);
    R("},");
    System.out.println("\tDone... " + (((double)(endTime - startTime))/1e9) + "\n");

    System.out.println("BFS...");
    {
      startTime = System.nanoTime();

      HashMap<Node, Boolean> found = new HashMap<Node, Boolean>((int)nv);
      HashMap<Node, Integer> distance = new HashMap<Node, Integer>((int)nv);
      for(int v = 1; v < nv; v++) {
	found.put(nodes[v], false);
	distance.put(nodes[v], -1);
      }
      found.put(nodes[0], true);
      distance.put(nodes[0], 0);

      LinkedList<Node> q = new LinkedList<Node>();
      q.add(nodes[0]);
      
      while(!q.isEmpty()) {
	Node v = q.poll();
	for(Relationship rel : v.getRelationships()) {
	  Node endNode = rel.getEndNode();
	  if(!found.get(endNode)) {
	    found.put(endNode, true);
	    distance.put(endNode, distance.get(v) + 1);
	    q.offer(endNode);
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
    double sssp_time = (((double)(endTime - startTime))/1e9);
    R("\"sssp\": {");
    R("\"name\":\"neo4j-std\",");
    R("\"time\":" +  sssp_time);
    R("},");
    System.out.println("\tDone... " + (((double)(endTime - startTime))/1e9) + "\n");

    System.out.println("Page Rank...");
    {
      double epsilon = 1e-8;
      double dampingfactor = 0.85;
      int maxiter = 100;

      startTime = System.nanoTime();
      HashMap<Node, Double> pr = new HashMap<Node, Double>((int)nv);
      HashMap<Node, Double> prTmp = new HashMap<Node, Double>((int)nv);

      for(Node v : GlobalGraphOperations.at(graphDb).getAllNodes()) {
	pr.put(v, 1.0/((double)nv));
      }

      int iter = maxiter;
      double delta = 1;

      while(delta > epsilon && iter > 0) {
	for(Node v : GlobalGraphOperations.at(graphDb).getAllNodes()) {
	  double myPrTmp = 0;
	  for(Relationship rel : v.getRelationships()) {
	    Node endNode = rel.getEndNode();
	    myPrTmp += pr.get(endNode) / ((double)countDegree(endNode));
	  }
	  prTmp.put(v, myPrTmp);
	}

	for(Node v : GlobalGraphOperations.at(graphDb).getAllNodes()) {
	  prTmp.put(v, prTmp.get(v) * dampingfactor + ((1.0-dampingfactor) / ((double)nv)));
	}

	delta = 0;
	for(Node v : GlobalGraphOperations.at(graphDb).getAllNodes()) {
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
    double pr_time = (((double)(endTime - startTime))/1e9);
    R("\"pr\": {");
    R("\"name\":\"neo4j-std\",");
    R("\"time\":" +  pr_time);
    R("},");
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

      Transaction tx = graphDb.beginTx();
      try {
	for(int a = 0; a < na; a++) {
	  long i = actions[2*a];
	  long j = actions[2*a+1];
	  if(i < 0) {
	    i =	~i;
	    j =	~j;
	    for(Relationship rel : nodes[(int)i].getRelationships()) {
	      if(nodes[(int)j] == rel.getEndNode()) {
		rel.delete();
		break;
	      }
	    }
	    for(Relationship rel : nodes[(int)j].getRelationships()) {
	      if(nodes[(int)i] == rel.getEndNode()) {
		rel.delete();
		break;
	      }
	    }
	  } else {
	    boolean handled = false;
	    for(Relationship rel : nodes[(int)i].getRelationships()) {
	      if(nodes[(int)j] == rel.getEndNode()) {
		rel.setProperty("weight", (Integer)rel.getProperty("weight") + 1);
		handled = true;
		break;
	      }
	    }
	    if(!handled) {
	      Relationship rel = nodes[(int)i].createRelationshipTo(nodes[(int)j], RelTypes.EDGE);
	      rel.setProperty("weight", 1);
	    }
	    handled = false;
	    for(Relationship rel : nodes[(int)j].getRelationships()) {
	      if(nodes[(int)i] == rel.getEndNode()) {
		rel.setProperty("weight", (Integer)rel.getProperty("weight") + 1);
		handled = true;
		break;
	      }
	    }
	    if(!handled) {
	      Relationship rel = nodes[(int)j].createRelationshipTo(nodes[(int)i], RelTypes.EDGE);
	      rel.setProperty("weight", 1);
	    }
	  }
	}

	tx.success();
      } catch (Exception e) {
	tx.failure();
      } finally {
	tx.finish();
      }

      endTime = System.nanoTime();
    }
    double eps = na/(((double)(endTime - startTime))/1e9);
    R("\"update\": {");
    R("\"name\":\"neo4j-std\",");
    R("\"time\":" +  eps);
    R("}");
    R("},");
    MemoryMXBean mem = ManagementFactory.getMemoryMXBean();
    long memory = (mem.getHeapMemoryUsage().getUsed() + mem.getNonHeapMemoryUsage().getUsed()) / 1024;

    R("\"na\":" + na + ",");
    R("\"mem\":" + memory);
    R("}");
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

  public static int countDegree(Node v) {
    int count = 0;
    for(Relationship r : v.getRelationships()) {
      count++;
    }
    return count;
  }

  private static void registerShutdownHook( final GraphDatabaseService graphDb ) {
    Runtime.getRuntime().addShutdownHook( 
      new Thread() {
	@Override
	public void run() {
	  graphDb.shutdown();
	}
      }
    );
  }
}
