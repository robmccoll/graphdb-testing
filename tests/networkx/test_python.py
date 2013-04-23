import sys
import struct
import networkx
import time
import json
import resource

# Reads STINGER binary format graph into networkx
# Uses pure Python functions to calculate connected components, SSSP, and PageRank
# Performs insertions and removals of all actions


class Timer:
  """http://preshing.com/20110924/timing-your-code-using-pythons-with-statement"""
  def start(self):
    self.start = time.time()
    return self

  def stop(self):
    self.end = time.time()
    self.interval = self.end - self.start


graphfile = open(sys.argv[1])
actionfile = open(sys.argv[2])

print 'reading graph and actions...'
t = Timer().start()

graphbin = graphfile.read()
actionbin = actionfile.read()

nv = struct.unpack('q',graphbin[8:16])[0];
ne = struct.unpack('q',graphbin[16:24])[0];
numactions = struct.unpack('q', actionbin[8:16])[0]

t.stop()

results = dict()
results["nv"] = nv;
results["ne"] = ne;
results["na"] = numactions;
results["type"] = "networkx";
results["results"] = dict();

print '   done ' + str(t.interval)

print '   nv : ' + str(nv)
print '   ne : ' + str(ne)
print '   na : ' + str(numactions)

print 'splitting csr...'

t = Timer().start()

off = struct.unpack(str(nv+1) + 'q', graphbin[24            :24+(nv+1)*8])
ind = struct.unpack(str(ne)   + 'q', graphbin[24+(nv+1)*8   :24+(nv+1+ne)*8])
wgt = struct.unpack(str(ne)   + 'q', graphbin[24+(nv+1+ne)*8:24+(nv+1+ne*2)*8])

t.stop()

results["results"]["build"] = dict();
results["results"]["build"]["name"] = "networkx-std";
results["results"]["build"]["time"] = t.interval
print '   done ' + str(t.interval)

print 'splitting actions...'

t = Timer().start()

actions = struct.unpack(str(numactions*2) + 'q', actionbin[16:16+numactions*2*8])

t.stop()

print '   done ' + str(t.interval)

print 'loading into networkx...'

t = Timer().start()

G = networkx.Graph()

for v in xrange(nv):
  for e in xrange(off[v],off[v+1]):
    G.add_edge(v,ind[e])

t.stop()

results["results"]["build"]["time"] += t.interval
print '   done ' + str(t.interval)

print 'connected components...'

t = Timer().start()

components = networkx.components.number_connected_components(G)

t.stop()

results["results"]["sv"] = dict();
results["results"]["sv"]["name"] = "networkx-std";
results["results"]["sv"]["time"] = t.interval

print '   done ' + str(t.interval)
print '   components ' + str(components)

print 'SSSP...'

t = Timer().start()

networkx.shortest_paths.unweighted.single_source_shortest_path_length(G, 0)

t.stop()

results["results"]["sssp"] = dict();
results["results"]["sssp"]["name"] = "networkx-std";
results["results"]["sssp"]["time"] = t.interval

print '   done ' + str(t.interval)

print 'PageRank...'

t = Timer().start()

prs = networkx.pagerank_alg.pagerank(G)

t.stop()

results["results"]["pr"] = dict();
results["results"]["pr"]["name"] = "networkx-std";
results["results"]["pr"]["time"] = t.interval

print '   done ' + str(t.interval)

print 'Insert / remove...'

t = Timer().start()

for a in xrange(numactions):
  if(actions[a*2] >= 0):
    G.add_edge(actions[a*2], actions[a*2+1])
  else:
    if(G.has_edge(~actions[a*2], ~actions[a*2+1])):
      G.remove_edge(~actions[a*2], ~actions[a*2+1])

t.stop()

results["results"]["update"] = dict();
results["results"]["update"]["name"] = "networkx-std";
results["results"]["update"]["time"] = t.interval

print '   done ' + str(t.interval)

results["mem"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss;

results_json = json.dumps(results)

for line in results_json.split('\n'):
  print("RSLT: " + line)
