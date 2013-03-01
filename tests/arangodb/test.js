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

  delta = 1;
  epsilon = 1e-8;
  dampingfactor = 0.85;
  maxiter = 100;
  iter = maxiter;
  tmp_pr = [];

  while(delta > epsilon && iter > 0) {
    for(i = 1; i <= nv; i++) {
      q = db.vertices.byExample({name : i});
      if(q.hasNext()) {
	v = q.next();
	tmp_pr[i] = 0;

	eq = db.edges.outEdges(v._id);
	if(eq) {
	  for(e in eq) {
	    edge = db.vertices.document(eq[e]._to);
	    tmp_pr[v.name] += (edge.pr / (db.edges.outEdges(eq[e]._to).length));
	  }
	}
      }
    }

    delta = 0;
    for(i = 1; i <= nv; i++) {
      q = db.vertices.byExample({name : i});
      if(q.hasNext()) {
	v = q.next();
	tmp_pr[v.name] = tmp_pr[v.name] * dampingfactor + ((1-dampingfactor) / (nv));

	my_delta = tmp_pr[v.name] - v.pr;
	delta += ((my_delta > 0) ? my_delta : (-my_delta));
	db.vertices.update(v, { "pr" : tmp_pr[v.name] });
      }
    }

    print("Iteration " + (maxiter - iter) + " delta " + delta);
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
page_rank(db, db.vertices.count());
var stop = new Date();
print("\tDone " + ((stop - start) / 1000));
