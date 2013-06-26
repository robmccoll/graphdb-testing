package spark.bagel.examples

import spark._
import spark.SparkContext._

import spark.bagel._
import spark.bagel.Bagel._

import scala.math._
import scala.collection.mutable.HashMap
import java.io._
import java.lang.management._

object App {

  class mVertex (
    var id : Int, 
    val componentId : Int, 
    var found : Boolean,
    var distance : Int,
    val pr : Double, 
    val prDelta : Double,
    val outEdges : HashMap[Int, Int], 
    val active : Boolean) 
      extends Vertex with Serializable

  class SVMessage (
    val targetId : Int,
    val componentId : Int)
      extends Message[Int] with Serializable

  class SVCombiner extends Combiner[SVMessage, Int] with Serializable {
    def createCombiner(msg: SVMessage): Int =
      msg.componentId
    def mergeMsg(combiner: Int, msg: SVMessage): Int = 
      min(msg.componentId, combiner)
    def mergeCombiners(a: Int, b: Int): Int =
      min(a, b)
  }

  class SPMessage (
    val targetId: Int, 
    val distance: Int)
      extends Message[Int] with Serializable

  class SPCombiner extends Combiner[SPMessage, Int] with Serializable {
    def createCombiner(msg: SPMessage): Int =
      msg.distance
    def mergeMsg(combiner: Int, msg: SPMessage): Int = 
      min(msg.distance, combiner)
    def mergeCombiners(a: Int, b: Int): Int =
      min(a, b)
  }

  class PRMessage (
    val targetId : Int,
    val pr: Double)
      extends Message[Int] with Serializable

  class PRAggregator extends spark.bagel.Aggregator[mVertex, Double] with Serializable {
      def createAggregator(vert: mVertex) : Double =
	vert.prDelta
      def mergeAggregators(a: Double, b: Double) : Double =
	a + b
  }

  class PRCombiner extends Combiner[PRMessage, Double] with Serializable {
    def createCombiner(msg: PRMessage): Double =
      msg.pr
    def mergeMsg(combiner: Double, msg: PRMessage): Double = 
      msg.pr + combiner
    def mergeCombiners(a: Double, b: Double): Double =
      a + b
  }

  class EdgeMessage (val targetId: Int, val destId: Int) extends Message[Int] with Serializable

  def R(result : String) {
    println("RSLT: " + result);
  }
  
  def main(args : Array[String]) {
    println( "Up and running..." )

    if(args.length < 2) {
      printf("Usage: make ARGS=\"graph.file actions.file\"\n")
      exit(-1)
    }

    printf("Reading graph from disk... %s\n", args(0))

    val graphIn = new FileInputStream(args(0))
    var startTime : Long = 0
    var endTime : Long = 0

    val endian = readLittleLong(graphIn)
    if(endian != 0x1234ABCD) {
      System.err.printf("Endianness failure " + endian + " \n")
      System.exit(-1)
    }

    val nv = readLittleLong(graphIn)
    val ne = readLittleLong(graphIn)
    val parts = 1

    print("\tGraph size nv " + nv + " ne " + ne)

    val off = (0 until nv+1).map(v => readLittleLong(graphIn))
    val ind = (0 until ne).map(v => readLittleLong(graphIn))
    val wgt = (0 until ne).map(v => readLittleLong(graphIn))

    println("\tDone\n")

    R("{");
    R("\"type\":\"bagel\",");
    R("\"nv\":" + nv + ",");
    R("\"ne\":" + ne + ",");
    R("\"results\": {");

    val host = args(2)
    //val sc = new SparkContext(host, "bageltest", System.getenv("SPARK_HOME"), Seq(System.getenv("SPARK_EXAMPLES_JAR")), sys.env)
    val sc = new SparkContext(host, "bageltest")

    startTime = System.nanoTime();

    /* holy syntactic sugar batman! - convert CSR into array of vertices with edge sequences*/
    //def toMap(start : Int, stop : Int) : HashMap[Int, Int] = { val 
    val verts = sc.parallelize((0 until nv).map{ 
      v : Int => (v, new mVertex(
	v, v, v == 0, if(v == 0) 0 else Int.MaxValue, 0, 100, { val edges = new HashMap[Int, Int](); (off(v) until off(v+1)).map{e => edges.put(ind(e), wgt(e))}; edges}, true))})

    endTime = System.nanoTime();

    var build_time = (((endTime - startTime).asInstanceOf[Double])/1e9);
    R("\"build\": {");
    R("\"name\":\"bagel\",");
    R("\"time\":" +  build_time);
    R("},");

    startTime = System.nanoTime();
    val emptySVMsgs = sc.parallelize(List[(Int, SVMessage)]())
    val resultSV = Bagel.run(sc, verts, emptySVMsgs, parts)(shiloachVishkin)
    val resultSVWithCombiner = Bagel.run(sc, verts, emptySVMsgs, combiner = new SVCombiner(), numPartitions = parts)(shiloachVishkinWithCombiner)
    endTime = System.nanoTime();
    var sv_time = (((endTime - startTime).asInstanceOf[Double])/1e9);
    R("\"sv\": {");
    R("\"name\":\"bagel\",");
    R("\"time\":" +  sv_time);
    R("},");

    startTime = System.nanoTime();
    val emptySPMsgs = sc.parallelize(List[(Int, SPMessage)]())
    val resultSP = Bagel.run(sc, resultSV, emptySPMsgs, parts)(shortestPaths)
    val resultSPWithCombiner = Bagel.run(sc, resultSV, emptySPMsgs, combiner = new SPCombiner(), numPartitions = parts)(shortestPathsWithCombiner)
    endTime = System.nanoTime();
    var sssp_time = (((endTime - startTime).asInstanceOf[Double])/1e9);
    R("\"sssp\": {");
    R("\"name\":\"bagel\",");
    R("\"time\":" +  sssp_time);
    R("},");

    startTime = System.nanoTime();
    val emptyPRMsgs = sc.parallelize(List[(Int, PRMessage)]())
    val resultPR = Bagel.run(sc, resultSV, emptyPRMsgs, combiner = new DefaultCombiner[PRMessage](), aggregator = Option(new PRAggregator()), numPartitions = parts, partitioner = new HashPartitioner(parts))(pageRank(1e-8, 0.85, 100, nv))
    val resultPRWithCombiner = Bagel.run(sc, resultSV, emptyPRMsgs, combiner = new PRCombiner(), aggregator = Option(new PRAggregator()), numPartitions = parts, partitioner = new HashPartitioner(parts))(pageRankWithCombiner(1e-8, 0.85, 100, nv))
    endTime = System.nanoTime();
    var pr_time = (((endTime - startTime).asInstanceOf[Double])/1e9);
    R("\"pr\": {");
    R("\"name\":\"bagel\",");
    R("\"time\":" +  pr_time);
    R("},");

    val actionIn = new FileInputStream(args(1))

    val actionendian = readLittleLong(actionIn)
    if(actionendian != 0x1234ABCD) {
      System.err.printf("Endianness failure " + actionendian + " \n")
      System.exit(-1)
    }

    val na = readLittleLong(actionIn)


    println("\t na " + na)

    val edgeActions = new Array[(Int, EdgeMessage)](na*2)
    (0 until na).foreach { a => { 
	val i = readLittleLong(actionIn)
	val j = readLittleLong(actionIn)
	edgeActions(2*a)   = (babs(i), new EdgeMessage(babs(i),j));
	edgeActions(2*a+1) = (babs(j), new EdgeMessage(babs(j),i));
    }}

    startTime = System.nanoTime();
    val resultInsertRemove = Bagel.run(sc, verts, sc.parallelize(edgeActions), parts)(insertRemoveEdges)
    endTime = System.nanoTime();
    var eps = na.asInstanceOf[Double] / (((endTime - startTime).asInstanceOf[Double])/1e9);
    R("\"update\": {");
    R("\"name\":\"bagel\",");
    R("\"time\":" +  eps);
    R("}");
    R("},");
    val mem = ManagementFactory.getMemoryMXBean();
    val memory = (mem.getHeapMemoryUsage().getUsed() + mem.getNonHeapMemoryUsage().getUsed()) / 1024;

    R("\"na\":" + na + ",");
    R("\"mem\":" + memory);
    R("}");
  }

  def shiloachVishkin(self : mVertex, msgs : Option[Array[SVMessage]], step : Int) : (mVertex, Array[SVMessage]) = {
    val newCid = min(msgs.getOrElse(Array(new SVMessage(self.id,self.componentId))).map(_.componentId).min, self.componentId)
    val halt = (newCid == self.componentId && step != 0)
    val msgout = if(!halt) self.outEdges.keySet.map(e => new SVMessage(e, newCid)).toArray else Array[SVMessage]()
    (new mVertex(self.id, newCid, self.found, self.distance, self.pr, self.prDelta, self.outEdges, !halt), msgout)
  }
  
  def shiloachVishkinWithCombiner(self : mVertex, msgmin : Option[Int], step : Int) : (mVertex, Array[SVMessage]) = {
    val newCid = min(msgmin.getOrElse(self.componentId), self.componentId)
    val halt = (newCid == self.componentId && step != 0)
    val msgout = if(!halt) self.outEdges.keySet.map(e => new SVMessage(e, newCid)).toArray else Array[SVMessage]()
    (new mVertex(self.id, newCid, self.found, self.distance, self.pr, self.prDelta, self.outEdges, !halt), msgout)
  }

  def shortestPaths(self : mVertex, msgin : Option[Array[SPMessage]], step : Int) : (mVertex, Array[SPMessage]) = {
    if(msgin.isDefined && !self.found) {
      (new mVertex(self.id, self.componentId, true, msgin.get(0).distance, self.pr, self.prDelta, self.outEdges, false), self.outEdges.keySet.map{ e => new SPMessage(e, msgin.get(0).distance+1)}.toArray)
    } else if(step == 0 && self.found) {
      (new mVertex(self.id, self.componentId, true, 0, self.pr, self.prDelta, self.outEdges, false), self.outEdges.keySet.map{ e => new SPMessage(e, 1)}.toArray)
    } else {
      (new mVertex(self.id, self.componentId, self.found, self.distance, self.pr, self.prDelta, self.outEdges, !self.found && self.componentId == 0), Array[SPMessage]())
    }
  }

  def shortestPathsWithCombiner(self : mVertex, msgin : Option[Int], step : Int) : (mVertex, Array[SPMessage]) = {
    if(msgin.isDefined && !self.found) {
      (new mVertex(self.id, self.componentId, true, msgin.get, self.pr, self.prDelta, self.outEdges, false), self.outEdges.keySet.map{ e => new SPMessage(e, msgin.get+1)}.toArray)
    } else if(step == 0 && self.found) {
      (new mVertex(self.id, self.componentId, true, 0, self.pr, self.prDelta, self.outEdges, false), self.outEdges.keySet.map{ e => new SPMessage(e, 1)}.toArray)
    } else {
      (new mVertex(self.id, self.componentId, self.found, self.distance, self.pr, self.prDelta, self.outEdges, !self.found && self.componentId == 0), Array[SPMessage]())
    }
  }

  def pageRank(epsilon : Double, damping : Double, maxIter : Int, nv : Double)(self: mVertex, msgin : Option[Array[PRMessage]], err : Option[Double], step : Int) :(mVertex, Array[PRMessage]) = {
    val prNew = if(msgin.isDefined) ((1 - damping) / nv + msgin.get.map(_.pr).sum * damping) else (1.0/nv)
    if(step < maxIter && (!err.isDefined || err.get > epsilon))
      (new mVertex(self.id, self.componentId, self.found, self.distance, prNew, abs(self.pr - prNew), self.outEdges, true), self.outEdges.keySet.map( e => new PRMessage(e, prNew / self.outEdges.size)).toArray)
    else
      (new mVertex(self.id, self.componentId, self.found, self.distance, prNew, abs(self.pr - prNew), self.outEdges, false), Array[PRMessage]())
  }

  def pageRankWithCombiner(epsilon : Double, damping : Double, maxIter : Int, nv : Double)(self: mVertex, msgin : Option[Double], err : Option[Double], step : Int) :(mVertex, Array[PRMessage]) = {
    val prNew = if(msgin.isDefined) ((1 - damping) / nv + msgin.get * damping) else (1.0/nv)
    if(step < maxIter && (!err.isDefined || err.get < epsilon))
      (new mVertex(self.id, self.componentId, self.found, self.distance, prNew, abs(self.pr - prNew), self.outEdges, true), self.outEdges.keySet.map( e => new PRMessage(e, prNew / self.outEdges.size)).toArray)
    else
      (new mVertex(self.id, self.componentId, self.found, self.distance, prNew, abs(self.pr - prNew), self.outEdges, false), Array[PRMessage]())
  }

  def insertRemoveEdges(self : mVertex, msgs : Option[Array[EdgeMessage]], step : Int) : (mVertex, Array[EdgeMessage]) = {
    if(msgs.isDefined)
      msgs.get.foreach( m => { 
	if(m.destId < 0) {
	  self.outEdges.remove(babs(m.destId))
	} else {
	  self.outEdges.put(m.destId, 1 + self.outEdges.getOrElse(m.destId, 1))
	}
      })
    (new mVertex(self.id, self.componentId, self.found, self.distance, self.pr, self.prDelta, self.outEdges, false), Array[EdgeMessage]())
  }

  def readLittleLong(in : FileInputStream) : Int = {
    val b = new Array[Byte](8)
    in.read(b)
    return (((b(3).asInstanceOf[Int] & 255) << 24) +
            ((b(2).asInstanceOf[Int] & 255) << 16) +
            ((b(1).asInstanceOf[Int] & 255) << 8 ) +
            ((b(0).asInstanceOf[Int] & 255) << 0 ))
  }

  def babs( i : Int) : Int = { if(i < 0) ~i else i }

}
