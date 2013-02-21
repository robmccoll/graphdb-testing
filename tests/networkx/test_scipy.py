import sys
import struct
import networkx
import time
import scipy

# Reads STINGER binary format graph into networkx
# Uses pure Python functions to calculate connected components, SSSP, and PageRank
# Performs insertions and removals of all actions

class Timer:
  """http://preshing.com/20110924/timing-your-code-using-pythons-with-statement"""
  def start(self):
    self.start = time.clock()
    return self

  def stop(self):
    self.end = time.clock()
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

print '   done ' + str(t.interval)

print 'splitting actions...'

t = Timer().start()

actions = struct.unpack(str(numactions*2) + 'q', actionbin[16:16+numactions*2*8])

t.stop()

print '   done ' + str(t.interval)

print 'loading into networkx...'

t = Timer().start()

A = scipy.sparse.csr_matrix((wgt,ind,off))
G = networkx.from_scipy_sparse_matrix(A)

t.stop()

print '   done ' + str(t.interval)

print 'connected components...'

t = Timer().start()

components = networkx.components.number_connected_components(G)

t.stop()

print '   done ' + str(t.interval)
print '   components ' + str(components)

print 'SSSP...'

t = Timer().start()

components = networkx.shortest_paths.unweighted.single_source_shortest_path_length(G, 0)

t.stop()

print '   done ' + str(t.interval)

print 'PageRank...'

t = Timer().start()

prs = networkx.pagerank_alg.pagerank_scipy(G)

t.stop()

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

print '   done ' + str(t.interval)

done = raw_input()
