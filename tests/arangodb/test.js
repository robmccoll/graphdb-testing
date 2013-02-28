function reset_values(db, nv) {
  q = db.vertices.all();

  while(q.hasNext()) {
    v = q.next();
    db.vertices.update(v, { "component" : v.name, "distance" : -1, "pr" : 1.0 / nv});
  }
}

function shiloach_vishkin(db) {
  changed = true;

  while(changed) {
    changed = false;
    eq = db.edges.all();

    while(eq.hasNext()) {
      e = eq.next();
      src = db._document(e._from);
      dst = db._document(e._to);
      if(src.component < dst.component) {
	db.vertices.update(src, { component : dst.component })
	changed = true;
      }
    }
  }
}

function count_components(db) {
  count = 0;
  q = db.vertices.all();

  while(q.hasNext()) {
    v = q.next();
    count += (v.name == v.component);
  }
  return count;
}

function bfs(db, start) {
  queue = [];
  src = db.vertices.byExample({name : start}).next();
  queue.push(src._id);

  while(queue.length) {
    v = queue.shift();
    v_doc = db.vertices.document(v)
    v_dist = v_doc.distance
    eq = db.edges.outEdges(v);
    for(e in eq) {
      edge = db.vertices.document(eq[e]._to);
      if(edge.distance == -1) {
	db.vertices.update(edge, { distance : v_dist + 1 })
	queue.push(edge._id)
      }
    }
  }
}

function page_rank(db, nv) {

  q = db.vertices.all();
  while(q.hasNext()) {
    v = q.next();
    db.vertices.update(v, { "component" : v.name, "distance" : -1, "pr" : 1.0 / nv, "pr_tmp" : 0});
  }

  delta = 1;
  epsilon = 1e-8;
  dampingfactor = 0.85;
  maxiter = 100;
  tmp_pr = [];

  while(delta > epsilon && iter > 0) {
    q = db.vertices.all();
    while(q.hasNext()) {
      v = q.next();
      tmp_pr[v.name] = 0;

      eq = db.edges.outEdges(v._id);
      for(e in eq) {
	edge = db.vertices.document(eq[e]._to);
	// TODO YOU ARE HERE
      }
    }

/*
    for(; vtxIt != vtxEnd; ++vtxIt++) {
      tmp_pr[*vtxIt] = 0;

      Graph::out_edge_iterator edgesIt, edgesEnd; 

      tie(edgesIt, edgesEnd) = out_edges(Graph::vertex_descriptor(*vtxIt), g);
      for(; edgesIt != edgesEnd; ++edgesIt) {
	tmp_pr[source(*edgesIt,g)] += (((double)pr[target(*edgesIt, g)]) / 
	  ((double) out_degree(target(*edgesIt, g), g)));
      }
    }

    for(uint64_t v = 0; v < nv; v++) {
      tmp_pr[v] = tmp_pr[v] * dampingfactor + (((double)(1-dampingfactor)) / ((double)nv));
    }

    delta = 0;
    for(uint64_t v = 0; v < nv; v++) {
      double mydelta = tmp_pr[v] - pr[v];

      if(mydelta < 0)
	mydelta = -mydelta;

      delta += mydelta;
      pr[v] = tmp_pr[v];
    }
    */

    iter--;
  }
}


print("Resetting values...");
reset_values(db, db.vertices.count());

print("Starting Shiloach-Vishkin...");
var start = new Date();
shiloach_vishkin(db);
count = count_components(db);
var stop = new Date();
print("\tDone " + ((stop - start) / 1000) + " Components " + count);

print("BFS...");
var start = new Date();
bfs(db, 1);
var stop = new Date();
print("\tDone " + ((stop - start) / 1000));

print("PageRank...");
var start = new Date();
page_rank(db, 1);
var stop = new Date();
print("\tDone " + ((stop - start) / 1000));
