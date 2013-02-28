<?php

namespace triagens\ArangoDb;

require 'init.php';

function get_next_int64($fp) {
  return unpack("L", fread($fp, 4))[1] + unpack("L", fread($fp, 4))[1] * 0x100000000;
}

function get_int64_array($fp, $size) {
  $return = array();
  for($i = 0; $i < $size; $i++) {
    $return[$i] = get_next_int64($fp);
  }
  return $return;
}

function new_vertex($documents, $id, $nv) {
  $v = new Document();
  $v->set("name", $id + 1);
  $v->set("component", $id);
  $v->set("distance", -1);
  $v->set("pr", 1.0/$nv);
  $documents->add("vertices",$v);
  return $v->getHandle();
}

if(count($argv) < 3) {
  $prog = $argv[0];
  print("Usage $prog graphfile actionsfile");
  exit();
}

print("Setting up DB connection...\n");

$connection	= new Connection($connectionOptions);
$collections	= new CollectionHandler($connection);
$documents	= new DocumentHandler($connection);
$edges		= new EdgeHandler($connection);

/* drop vertex if it exists */
try {
  $vertex = $collections->get("vertices");
  $collections->drop($vertex);
} catch(Exception $e) {
  // do nothing
}

/* drop edge if it exists */
try { 
  $edge = $collections->get("edges");
  $collections->drop($edge);
} catch(Exception $e) {
  // do nothing
}

$vertex = new Collection();
$vertex->setName("vertices");
$collections->add($vertex);
$vertex = $vertex->getId();

$edge = new Collection();
$edge->setName("edges");
$edge->setType(3);
$collections->add($edge);
$edge = $edge->getId();

print("Reading graph file...\n");

$graphfilename = $argv[1];
$actionsfilename = $argv[2];

$fp = fopen($graphfilename, "rb");

if(0x1234ABCD != get_next_int64($fp)) {
  print("Endianness check failed.  Endianness swap not implemented.");
  exit();
}

$nv = get_next_int64($fp);
$ne = get_next_int64($fp);

$off = get_int64_array($fp, $nv+1);
$ind = get_int64_array($fp, $ne);
$wgt = get_int64_array($fp, $ne);

fclose($fp);

print("Inserting into DB...\n");

//$batch = new Batch($connection);
//$batch->activate();

$vertices = array_fill(0, $nv, "-1");
$time_start = microtime(true);

for($v = 0; $v < $nv; $v++) {
  $s = $vertices[$v];
  if($s == "-1") {
    $s = new_vertex($documents, $u, $nv);
    $vertices[$v] = $s;
  }
  for($i = $off[$v]; $i < $off[$v+1]; $i++) {
    $u = $ind[$i];
    $d = $vertices[$u];
    if($d == "-1") {
      $d = new_vertex($documents, $u, $nv);
      $vertices[$u] = $d;
    }
    $edg = new Edge();
    $edg->set("weight", $wgt[$i]);
    $edges->saveEdge($edge, $s, $d, $edg, true);
  }
}

//$batch->process();

?>
