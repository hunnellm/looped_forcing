from collections import deque

def gzerosgame(g, F=None, B=None):
	"""
	Return the derived set for a given graph g with set of banned edges B and a initial set of vertices. The derived set is given by doing generalized zero forcing process. That is, if y is the only white neighbor of x and xy is not banned, then x could force y into black.

	Input:
		g: a simple graph
		F: a list of vertices of g
		B: a list of tuples representing banned edges of g

	Output:
		A set of black vertices when zero forcing process stops.

	Examples:
		sage: gzerosgame(graphs.PathGraph(5),[0])
		set([0, 1, 2, 3, 4])
		sage: gzerosgame(graphs.PathGraph(5),[0],[(1,2)])
		set([0, 1])
	"""
	if F is None:
		F = []
	if B is None:
		B = []
	# Precompute neighbor sets once for O(1) lookup
	neighbors = {v: set(g.neighbors(v)) for v in g.vertices()}
	# Normalize banned edges to frozensets for O(1) membership, handling both directions
	banned = {frozenset(e) for e in B}
	black = set(F)
	in_queue = set(black)
	queue = deque(black)
	while queue:
		x = queue.popleft()
		in_queue.discard(x)
		whites = neighbors[x] - black
		if len(whites) != 1:
			continue
		(y,) = whites
		if frozenset({x, y}) in banned:
			continue
		# x forces y
		black.add(y)
		# y is now black: add y and all black neighbors of y back to the worklist
		# because their white-neighbor counts may have changed
		if y not in in_queue:
			queue.append(y)
			in_queue.add(y)
		for z in neighbors[y]:
			if z in black and z not in in_queue:
				queue.append(z)
				in_queue.add(z)
	return black

def gZ_leq(graph, support=None, bannedset=None, i=None):
	"""
	For a given graph with support and banned set, if there is a zero forcing set of size i then return it; otherwise return False.

	Input:
		graph: a simple graph
		support: a list of vertices of g
		bannedset: a list of tuples representing banned edges of graph
		i: an integer, the function check gZ <= i or not

	Output:
		if F is a zero forcing set of size i and support is a subset of F, then return F
		False otherwise

	Examples:
		sage: gZ_leq(graphs.PathGraph(5),[],[],1)
		set([0])
		sage: gZ_leq(graphs.PathGraph(5),[],[(0,1)],1) 
		False
	"""
	if support is None:
		support = []
	if bannedset is None:
		bannedset = []
	if i < len(support):
		return False
	j = i - len(support) # additional number of black vertices
	support_set = set(support)
	# Use list comprehension to avoid O(n^2) list.remove calls
	VX = [v for v in graph.vertices() if v not in support_set]
	order = graph.order()
	for subset in Subsets(VX, j):
		test_set = support_set.union(subset) # the set is tested to be a zero forcing set
		outcome = gzerosgame(graph, test_set, bannedset)
		if len(outcome) == order:
			return test_set
	return False

def find_gzfs(graph, support=None, bannedset=None, upper_bound=None, lower_bound=None):
	"""
	For a given graph with support and banned set, return the an optimal generalized zero forcing set. If upper_bound is less than the generalized zero forcing number then return ['wrong']. If lower_bound is greater than the generalized zero forcing number then the return value will not be correct

	Input:
		graph: a simple graph
		support: a list of vertices of g
		bannedset: a list of tuples representing banned edges of graph
		upper_bound: an integer supposed to be an upper bound of gZ. 
		lower_bound: an integer supposed to be a lower bound of gZ. The two bounds may shorten the computation time. But one may leave it as default value if one is not sure.

	Output:
		if F is an optimal zero forcing set of size i then return F. If upper_bound is less than the general zero forcing number then return ['wrong'].

	Examples:
		sage: find_gzfs(graphs.PathGraph(5))
		set([0])
		sage: find_gzfs(graphs.PathGraph(5),[1],[(3,2)])
		set([0, 1, 3])
	"""
	if support is None:
		support = []
	if bannedset is None:
		bannedset = []
	support_set = set(support)
	# Use list comprehension to avoid O(n^2) list.remove calls
	VX = [v for v in graph.vertices() if v not in support_set]
	order = graph.order()
	s = len(support)
	if upper_bound is None:
		upper_bound = order # the default upper bound
	if lower_bound is None:
		lower_bound = len(VX) # temporary lower bound
		for v in VX:
			N = set(graph.neighbors(v))
			D = N.difference(support_set)
			lower_bound = min([lower_bound, len(D)])
		for v in support:
			N = set(graph.neighbors(v))
			D = N.difference(support_set)
			lower_bound = min([lower_bound, len(D) - 1])
		lower_bound = lower_bound + s # the default lower bound
	i = upper_bound
	find = 1 # does sage find a zero forcing set of size i
	outcome = ['wrong'] # default outcome
	while i >= lower_bound and find == 1:
		find = 0
		leq = gZ_leq(graph, support, bannedset, i) # check gZ <= i or not
		if leq != False:
			outcome = leq
			find = 1
			i = i - 1
	return outcome

def find_gZ(graph, support=None, bannedset=None, upper_bound=None, lower_bound=None):
	"""
	For a given graph with support and banned set, return the zero. upper_bound and lower_bound could be left as default value if one is not sure.

	Input:
		graph: a simple graph
		support: a list of vertices of g
		bannedset: a list of tuples representing banned edges of graph
		upper_bound: an integer supposed to be an upper bound of gZ. 
		lower_bound: an integer supposed to be a lower bound of gZ. The two bounds may shorten the computation time. But one may leave it as default value if one is not sure.

	Output:
		the generalized zero forcing number

	Examples:
		sage: find_gZ(graphs.PathGraph(5))            
		1
		sage: find_gZ(graphs.PathGraph(5),[1],[(3,2)])
		3
	"""
	return len(find_gzfs(graph, support, bannedset, upper_bound, lower_bound))

def X(g):
	"""
	For a given graph g, return the verices set X of a part of the bipartite used to compute the exhaustive zero forcing number.

	Input:
		g: a simple graph

	Output:
		a list of tuples ('a',i) for all vertices i of g

	Examples:
		sage: X(graphs.PathGraph(5))
		[('a', 0), ('a', 1), ('a', 2), ('a', 3), ('a', 4)]
	"""
	return [('a',i) for i in g.vertices()]

def Y(g):
	"""
	For a given graph g, return the verices set Y of the other part of the bipartite used to compute the exhaustive zero forcing number.

	Input:
		g: a simple graph

	Output:
		a list of tuples ('b',i) for all vertices i of g

	Examples:
		sage: Y(graphs.PathGraph(5))
		[('b', 0), ('b', 1), ('b', 2), ('b', 3), ('b', 4)]
	"""
	return [('b',i) for i in g.vertices()]

def _tilde_bipartite_base(g):
	"""
	Build the invariant (I-independent) portion of the tilde_bipartite graph for g.
	This includes all vertices X(g) and Y(g) and all adjacency-induced edges,
	but none of the bridge edges (('a',i),('b',i)) that depend on I.
	"""
	h = Graph()
	h.add_vertices(X(g))
	h.add_vertices(Y(g))
	for i in g.vertices():
		for j in g.neighbors(i):
			h.add_edge(('a', i), ('b', j))
	return h

def _tilde_bipartite_with_I(base, I):
	"""
	Given a base bipartite graph (from _tilde_bipartite_base) and an index set I,
	return the full tilde_bipartite graph by copying base and adding bridge edges.
	"""
	h = base.copy()
	for i in I:
		h.add_edge(('a', i), ('b', i))
	return h

def tilde_bipartite(g, I=None):
	"""
	For a given graph g and an index set I, return the bipartite graph \widetilde{G}_I used to compute the exhaustive zero forcing number.

	Input:
		g: a simple graph
		I: a list of vertices of g

	Output:
		the bipartite graph \widetilde{G}_I

	Examples:
		sage: h=tilde_bipartite(graphs.PathGraph(5),[1])
		sage: h.vertices()
		[('a', 0), ('a', 1), ('a', 2), ('a', 3), ('a', 4), ('b', 0), ('b', 1), ('b', 2), ('b', 3), ('b', 4)]
		sage: h.edges()
		[(('a', 0), ('b', 1), None), (('a', 1), ('b', 0), None), (('a', 1), ('b', 1), None), (('a', 1), ('b', 2), None), (('a', 2), ('b', 1), None), (('a', 2), ('b', 3), None), (('a', 3), ('b', 2), None), (('a', 3), ('b', 4), None), (('a', 4), ('b', 3), None)]
	"""
	if I is None:
		I = []
	return _tilde_bipartite_with_I(_tilde_bipartite_base(g), I)

def find_EZ(g,bound=None):
	"""
	For a given graph g, return the exhaustive zero forcing number of g. A given bound may shorten the computation.

	Input:
		g: a simple graph
		bound: a integer as an upper bound. It could be left as default value if one is not sure.

	Output:
		the exhaustive zero forcing number (EZ) of g

	Examples:
		sage: find_EZ(graphs.PathGraph(5))
		1
		sage: h=graphs.CycleGraph(5)
		sage: h.add_vertices([5,6,7,8,9])
		sage: h.add_edges([(0,5),(1,6),(2,7),(3,8),(4,9)])
		sage: find_EZ(h) # the EZ of a 5-sun
		2
	"""
	order=g.order()
	Z=find_gZ(g) # without support and banned set, the value is the original zero forcing number
	if bound==None:
		bound=Z # default upper bound
	gZ_bound=bound+order
	V=set(g.vertices())
	e=-1 # temporary output
	# Build the invariant base graph once; add bridge edges per subset I
	base = _tilde_bipartite_base(g)
	Yg = Y(g)
	for I in Subsets(V):
		h = _tilde_bipartite_with_I(base, list(I))
		leq=gZ_leq(h,Yg,[],e) # this avoid abundant computation
		if leq==False:
			e=find_gZ(h,Yg,[],gZ_bound,e+1)
			# in this case, we already know e+1-order<=gZ-order<=bound and so e+1<=gZ<=gZ_bound
		if e==gZ_bound:
			break
	return e-order # EZ=max-order
	
def bridged_edges(J):
    """
    For a give subset J of vertices, return the corresponding edges between X and Y in tilde_bipartite.
    
    Input:
        J: a subset of vertices.
        
    Output:
        [(("a",j),("b",j)) for j in J].
        
    Examples:
        sage: J=[1,3,5];
        sage: print bridged_edges(J);
        [(('a', 1), ('b', 1)), (('a', 3), ('b', 3)), (('a', 5), ('b', 5))]
    """        
    return [(("a",j),("b",j)) for j in J];

def find_loopedZ(g, I, J=None, _base=None):
    """
    For a given graph g and the index of the vertices with loops, return the zero forcing number of this looped graph.
    
    Input:
        g: a simple graph, the underlying graph of the looped graph.
        I: the index of the vertices with exactly one loop.
        J: the index of the vertices with double loops.
        _base: optional precomputed base bipartite graph (from _tilde_bipartite_base) for reuse.
    
    Output:
        the zero forcing number of this (multi-)looped graph.
        
    Examples:
        sage: g = Graph({0:[1],1:[0]});
        sage: I=[0,1];
        sage: find_loopedZ(g,I)
        1
        sage: g = Graph({0:[1],1:[0]});
        sage: I=[0];
        sage: find_loopedZ(g,I)
        0   
        sage: g = Graph({0:[1],1:[0]});
        sage: I=[0];
        sage: J=[1];
        sage: find_loopedZ(g,I,J) 
        1
    """
    if J is None:
        J = []
    if _base is None:
        _base = _tilde_bipartite_base(g)
    return find_gZ(_tilde_bipartite_with_I(_base, I + J), Y(g), bridged_edges(J)) - g.order()
    
def diagonal_analysis(g, Z=None):
    """
    For a given graph, do the diagonal analysis and return the set of zero-vertices and nonzero-vertices.
    Input:
        g: simple graph considered.
        Z: if no input, then Z=zero forcing number; occasionally it can also be replaced by other value.
        
    Output:
        A dictionary of vertices and its zero-nonzero pattern: 0=zero, 1=nonzero, 2=free, and -1=impossible.
        
    Examples:
        sage: diagonal_analysis(graphs.CompleteGraph(3));
        {0: 1, 1: 1, 2: 1}
        sage: g = Graph({0:[1],1:[0,2,5],2:[1,3,6],3:[2,4,7],4:[3,5,8],5:[1,4,9],6:[2],7:[3],8:[4],9:[5]});
        sage: diagonal_analysis(g);
        #This is the 5-sun.
        {0: -1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: -1, 7: -1, 8: -1, 9: -1}
        sage: g=graphs.CompleteBipartiteGraph(1,3);
        sage: diagonal_analysis(g);
        {0: 2, 1: 0, 2: 0, 3: 0}
    """
    if Z is None:
        Z = find_gZ(g)
    # Build the base bipartite graph once and reuse across all vertices
    base = _tilde_bipartite_base(g)
    Vg = list(g.vertices())
    diag = {}
    for v in Vg:
        both = 0
        diag[v] = 2
        J = [u for u in Vg if u != v]
        if find_loopedZ(g, [], J, _base=base) < Z:
            both += 1
            diag[v] = 1
        if find_loopedZ(g, [v], J, _base=base) < Z:
            both += 1
            diag[v] = 0
        if both == 2:
            diag[v] = -1
    return diag
