from collections import deque

def find_Z(g, support=None, ban=None, return_sets=False):
    """
    Compute the minimum size of a zero forcing set for CCR-Z (Z_game).

    Input:
        g: a simple graph
        support: optional list/set of vertices that must be included in every candidate set
        ban: optional list/set of banned vertices (passed to Z_game); banned vertices
             cannot make a force but can still be forced
        return_sets: bool (default False)
            * False – return the minimum size (an int)
            * True  – return a sorted list of lists, each inner list a minimum ZFS

    Output:
        int  when return_sets is False
        list when return_sets is True

    Examples:
        sage: g = graphs.PathGraph(5)
        sage: find_Z(g)
        1
        sage: find_Z(g, return_sets=True)
        [[0], [4]]
        sage: g = graphs.CompleteGraph(3)
        sage: find_Z(g)
        2
    """
    if support is None:
        support = []
    if ban is None:
        ban = []
    support_set = set(support)
    V = list(g.vertices())
    # Early impossible case
    if not support_set.issubset(set(V)):
        return [] if return_sets else False

    # Try sizes from |support| up to |V|
    for k in range(len(support_set), len(V) + 1):
        results = []
        for subset in Subsets([v for v in V if v not in support_set], k - len(support_set)):
            B = list(support_set.union(subset))
            derived = Z_game(g, B, ban=ban)
            if len(set(derived)) == g.order():
                if not return_sets:
                    return k
                results.append(sorted(B))
        if return_sets and results:
            # Sort deterministically
            results = sorted(results, key=lambda lst: tuple(lst))
            return results
    return [] if return_sets else False


def find_Zplus(g, support=None, return_sets=False):
    """
    Compute the minimum size of a zero forcing set for CCR-Zplus (Zplus_game).

    Input:
        g: a simple graph
        support: optional list/set of vertices that must be included in every candidate set
        return_sets: bool (default False)
            * False – return the minimum size (an int)
            * True  – return a sorted list of lists, each inner list a minimum Zplus ZFS

    Output:
        int  when return_sets is False
        list when return_sets is True

    Examples:
        sage: g = graphs.PathGraph(5)
        sage: find_Zplus(g)
        1
        sage: find_Zplus(g, return_sets=True)
        [[0], [4]]
        sage: g = graphs.CompleteGraph(4)
        sage: find_Zplus(g)
        3
    """
    if support is None:
        support = []
    support_set = set(support)
    V = list(g.vertices())
    if not support_set.issubset(set(V)):
        return [] if return_sets else False

    for k in range(len(support_set), len(V) + 1):
        results = []
        for subset in Subsets([v for v in V if v not in support_set], k - len(support_set)):
            B = list(support_set.union(subset))
            derived = Zplus_game(g, B)
            if len(set(derived)) == g.order():
                if not return_sets:
                    return k
                results.append(sorted(B))
        if return_sets and results:
            results = sorted(results, key=lambda lst: tuple(lst))
            return results
    return [] if return_sets else False


def find_Zell_zero(g, support=None, ban=None, return_sets=False):
    """
    Compute the minimum size of a zero forcing set for CCR-Zell (Zell_game) (i.e. no loops on isolates).

    Note: Zell_game internally deletes isolated vertices from a copy of g (since they
    cannot be affected by forces). A set B "forces all vertices" of the original g iff:
      - every isolated vertex of g is in B, and
      - Zell_game forces all non-isolated vertices.

    Input:
        g: a simple graph
        support: optional list/set of vertices that must be included in every candidate set
        ban: optional list/set of banned vertices (passed to Zell_game)
        return_sets: bool (default False)
            * False – return the minimum size (an int)
            * True  – return a sorted list of lists, each inner list a minimum Zell ZFS

    Output:
        int  when return_sets is False
        list when return_sets is True

    Examples:
        sage: g = graphs.PathGraph(5)
        sage: find_Zell_zero(g)
        1
        sage: find_Zell_zero(g, return_sets=True)
        [[0], [4]]
        sage: g = graphs.PathGraph(3)
        sage: g.add_vertex(3)  # isolated
        sage: find_Zell_zero(g)
        2
    """
    if support is None:
        support = []
    if ban is None:
        ban = []
    V = list(g.vertices())
    iso = {v for v in V if g.degree(v) == 0}

    # Must include all isolated vertices to "color all vertices" in original g
    support_set = set(support).union(iso)
    if not support_set.issubset(set(V)):
        return [] if return_sets else False

    # Build the non-isolated induced subgraph; Zell_game will do the same deletion,
    # but we need its order for the success criterion.
    non_iso = [v for v in V if v not in iso]
    if not non_iso:
        # Graph is entirely isolated vertices: the only way to color all is take all of them.
        allv = sorted(V)
        return [allv] if return_sets else len(V)

    h = g.subgraph(non_iso)
    target = h.order()

    for k in range(len(support_set), len(V) + 1):
        results = []
        for subset in Subsets([v for v in V if v not in support_set], k - len(support_set)):
            B = list(support_set.union(subset))
            derived = Zell_game(g, B, ban=ban)
            # derived may include isolated vertices if they were in B; remove them for the check
            forced_non_iso = set(derived).intersection(set(non_iso))
            if len(forced_non_iso) == target:
                if not return_sets:
                    return k
                results.append(sorted(B))
        if return_sets and results:
            results = sorted(results, key=lambda lst: tuple(lst))
            return results
    return [] if return_sets else False


def Z_game(g,B,ban=[]):
    """
    Input:
        g: a simple graph
        B: a set of initial blue vertices
        ban: a set of banned vertices
    Output:
        return the derived set under the regular CCR-Z.
        Note: vertices in ban cannot make a force, but are still white neighbors if white.
    """    
    V=g.vertices();
    white_neighbors={}; #a dictionary with the structure {v: list of white neighbors}
    white_numbers={}; #a dictionary with the structure {v: number of white neighbors}
    for v in V:
        nbh=g.neighbors(v);
        for b in B:
            try:
                nbh.remove(b);
            except ValueError:
                pass;
        white_neighbors[v]=nbh;
        white_numbers[v]=len(nbh);
    queue=copy(B); #queue stores list of vertices that can possibly make a force
    derived_set=copy(B); #derived_set stores the set of blue vertices
    whole_loop=True;
    while whole_loop: #keep searching if queue!=[]
        try:
            v=queue[0];
            queue.remove(v);
            if v not in ban and white_numbers[v]==1:
                u=white_neighbors[v][0]; #the only white neighbor
                derived_set.append(u); #make the force
                #update white_numbers, white_neighbors, and queue
                if white_numbers[u]==1: 
                    queue.append(u);
                u_nbr=g.neighbors(u);
                for w in u_nbr:
                    white_neighbors[w].remove(u);
                    white_numbers[w]+=-1;
                    if w in derived_set and white_numbers[w]==1:
                        queue.append(w);
        except IndexError:
            whole_loop=False;
    return derived_set;

def Zell_game(g,B,ban=[]): ##every non-isolated vertex can force itself if no other white neighbors.
    """
    Input:
        g: a simple graph
        B: a set of initial blue vertices
        ban: a set of banned vertices
    Output:
        return the derived set under the CCR-Zell.  
        That is, in addition to CCR-Z, every non-isolated vertex can force itself if 
        no other white neighbors.
        Note: vertices in ban cannot make a force, but are still white neighbors if white.
    """    
    #Build h from g by deleting all isolated vertices,
    # since they are not going to change the outcome. 
    h=g.copy();
    for v in g.vertices():
        if h.degree(v)==0:
            h.delete_vertex(v);
    V=h.vertices();
    white_neighbors={}; #a dictionary with the structure {v: list of white closed neighbors}
    white_numbers={}; #a dictionary with the structure {v: number of white closed neighbors}
    for v in V:
        nbh=h.neighbors(v);
        nbh.append(v);
        for b in B:
            try:
                nbh.remove(b);
            except ValueError:
                pass;
        white_neighbors[v]=nbh;
        white_numbers[v]=len(nbh);
    queue=copy(V); #queue stores list of vertices that can possibly make a force
    derived_set=copy(B); #derived_set stores the set of blue vertices
    whole_loop=True;
    while whole_loop: #keep searching if queue!=[]
        try:
            v=queue[0];
            queue.remove(v);
            if v not in ban and white_numbers[v]==1:
                u=white_neighbors[v][0]; #the only white neighbor
                derived_set.append(u); #make the force
                #update white_numbers, white_neighbors, and queue
                u_nbr=h.neighbors(u);
                u_nbr.append(u);
                for w in u_nbr:
                    white_neighbors[w].remove(u);
                    white_numbers[w]+=-1;
                    if white_numbers[w]==1:
                        queue.append(w);
        except IndexError:
            whole_loop=False;
    return derived_set;
    
def Zplus_game(g,B):
    """
    Input:
        g: a simple graph
        B: a set of initial blue vertices
    Output:
        return the derived set under the CCR-Zplus.  
        That is, apply CCR-Z to each white branch.
        Note: the banned set is not implemented, since Zvcoc has no Zplus version.
    """   
    white_graph=g.copy(); #recording the induced subgraph on white
    derived_set=copy(B); #derived_set stores the set of blue vertices
    for b in derived_set:
        white_graph.delete_vertex(b);
    whole_loop=True;
    while whole_loop: 
        whole_loop=False; #open again only when something found.
        whole_extra_B=[]; #extra blue vertices found in this round
        partition=white_graph.connected_components();
        for com in partition:
            considered_set=copy(com);
            for v in derived_set:
                considered_set.append(v);
            considered_graph=g.subgraph(considered_set); #this is a white branch
            extra_B=Z_game(considered_graph,derived_set); #apply Z_game to this white branch
            for v in derived_set:
                extra_B.remove(v);
            for v in extra_B:
                whole_loop=True; #something found, open whole_loop
                whole_extra_B.append(v);
        #update new derived_set and new white_graph.
        for v in whole_extra_B:
            derived_set.append(v);
            white_graph.delete_vertex(v);
    return derived_set;  
    


def pt_plus(G,S):
    """
    pt_plus(G,S) returns the propagation time for
    positive semidefinite zero forcing on the graph G
    given a set S

    
    INPUT:
            graph G, vertex subset S
        
    OUTPUT: 
            Returns the propagation time based on the
            positive semidefinite color change rule.
    
    EXAMPLES:
        sage: pt_plus(graphs.CompleteGraph(6),[0,1,2,3])
        -1
        sage: pt_plus(graphs.PetersenGraph(),[0,2,8,7]) 
        2
    """
    prop_time = 0
    H = copy(G)
    H.delete_vertices(S)
    T = list(copy(S))
    while len(H.connected_components())==1:
        prop_time += 1
        new_T=copy(T)
        for v in T:
            N = set(G.neighbors(v))
            N.difference_update(set(T))
            if len(N) == 1:
                new_T.append(N.pop())
        if T == new_T:
            return -1
        T = new_T
        H = copy(G)
        H.delete_vertices(T)
    if len(H.connected_components())==0:
        return prop_time
    else:
        sub_times=[]
        for V in H.connected_components():
            B=set(G.vertices())
            B.difference_update(set(V))
            B=list(B)
            sub_times.append(pt_plus(G,B))
            if sub_times[-1]==-1:
                return -1
        return prop_time+max(sub_times)
        
        

def pt_plus_set(G):
    """
    given a graph G computes all possible
    semidefinite propagation times and returns as
    list of lists
    
    INPUT:
            graph G
        
    OUTPUT: 
            A list of lists, where each sublist[i] is the set of propagation
            times for positive semidefinite zero forcing sets of size i.  If
            sublist[i] is empty then no positive semidefinite zero forcing set
            of size i exists.
    
    EXAMPLES:
        sage: pt_plus_set(graphs.CompleteGraph(4))
        [[],[],[1],[0]]
        sage: pt_plus_set(atlas_graphs[5]) #path on three vertices and an isolated vertex
        [[], [], [1, 2], [1], [0]]
    """
    H=copy(G)
    H.relabel()
    all_times=[[]]
    for i in range(1,H.order()):
        #the i^th entry of i_times corresponds to forcing sets of size i
        i_times=[]
        #checking the propagation time of every vertex subset of size i
        for S in combinations(range(H.order()),i):
            q=pt_plus(H,S)
            if q>-1 and q not in i_times:
                i_times.append(q)
        i_times.sort()
        all_times.append(i_times)
    all_times.append([0])
    return all_times
    

def Zplus_gen(G):
    """
    Given a graph G, connected or disconnected, calculates the
    positive semidefinite zero forcing number for each connected
    component and returns the sum.
    
    INPUT:
            graph G
        
    OUTPUT: 
            If the graph is connected returns the positive
    semidefinite zero forcing number.
            If the graph is disconneceted, returns the sum of the 
    positive semidefinite zero forcing number of each connected
    component.
    
    EXAMPLES:
        sage: Zplus_gen(graphs.CompleteGraph(6))
        5
        sage: Zplus_gen(atlas_graphs[5]) 
        2
    """
    
    z = 0
    #for each connected component find the positive semidefinite zero forcing number and sum
    for H in G.connected_components_subgraphs():
        #The positive semidefinite zero forcing number of a tree is 1 
        if H.is_tree():
            z = z+1
        else:
            z = z + Zplus(H)
        
    return z




def min_pt_plus_int(G):
    """
    Given a graph G this function computes the minimum and maximum propagation time
    of G, along with vertex subsets that realize them and the propagation time interval.
    
    INPUT:
            graph G
        
    OUTPUT: 
            Returns a list of three items: positive semidefinite propagation time interval, 
            a minimum positive semidefinite zero forcing set that realizes the minimum 
            propagation time, a minimum positive semidefinite zero forcing set that realizes
            the maximum propagation time. 
    
    EXAMPLES:
        sage: min_pt_plus_int(graphs.CompleteGraph(6))
        [[1],[0,1,2,3,4],[0,1,2,3,4]]
        sage: min_pt_plus_int(graphs.PetersenGraph()) 
        [[1, 2, 3], [0, 2, 8, 9], [0, 1, 2, 8]]
    """
    min_times = []
    max_set = []
    min_set = []
    H = copy(G)
    H.relabel()
    z=Zplus_gen(H)
    c = 0
    #checking the propagation time for all minimum positive semidefinite zero forcing sets.
    for S in Combinations(range(H.order()),z):
        q = pt_plus(H,S)
        if q > -1 and q not in min_times and c==0:
            c = c+1
            min_times.append(q)
            
        if q > -1 and c > 0:
            m = min(min_times)
            n = max(min_times)
            
            #adding sets that realize minimum prop time
            if q <= m:
                min_set.append(S)
                if q < m:
                    min_set = []
                    min_set.append(S) 
            #adding sets that realize maximum prop time
            if q >= n:
                max_set.append(S)
                if q > n:
                    max_set = []
                    max_set.append(S)
            if q not in min_times:
                min_times.append(q) 
        
    min_times.sort()
    return [min_times,min_set[0],max_set[0] ]
    
    

def full_pt_plus_int(G):
    """
        Given a graph G this function tells us whether or not
        the graph has a full positive semidefinite
        propagation time interval.
    
    INPUT:
            graph G
        
    OUTPUT: 
            Returns True if the propagation time interval is full, returns False else.
           
    
    EXAMPLES:
        sage: full_pt_plus_int(graphs.CompleteGraph(6))
        True
        sage: full_pt_plus_int(graphs.PetesonGraph) 
        True
    """
    z = Zplus_gen(G)
    #We know trees have a full propagation time interval
    if z == 1:
        return True
    else:
        p = min_pt_plus_int(G)[0] #the minimum prop time
        l = len(p) #how long the prop time interval is
        d = p[l-1] - p[0] # the difference between the minimum and maximum prop time
        D = (d + 1) - l 
    if D == 0:
        return True
    return False
    

def gzerosgame(g, F=None, B=None):
    """
    Return the derived set for a given graph g with set of banned edges B and a initial set of vertices. The derived set is given by doing generalized zero forcing process. That is, if y is the only whit[...] 

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

def gZ_leq_all(graph, support=None, bannedset=None, i=None):
    """
    For a given graph with support and banned set, return *all* zero forcing sets of size i
    (each a frozenset containing the support plus exactly i-len(support) additional vertices).
    Returns an empty list when no such set exists.

    Input:
        graph: a simple graph
        support: a list of vertices of graph
        bannedset: a list of tuples representing banned edges of graph
        i: an integer; collect every ZFS of exactly this size

    Output:
    A list of frozensets, each a zero forcing set of size i that includes the support.
    The list is in the iteration order produced by Subsets (stable/deterministic).

    Examples:
        sage: gZ_leq_all(graphs.PathGraph(5), [], [], 1)
        [frozenset({0})]
        sage: gZ_leq_all(graphs.PathGraph(5), [], [], 2)
        [frozenset({0, 1}), frozenset({0, 2}), frozenset({0, 3}), frozenset({0, 4}), frozenset({1, 2}), frozenset({1, 3}), frozenset({1, 4}), frozenset({2, 3}), frozenset({2, 4}), frozenset({3, 4})]
    """
    if support is None:
        support = []
    if bannedset is None:
        bannedset = []
    if i < len(support):
        return []
    j = i - len(support)  # additional vertices beyond support
    support_set = set(support)
    VX = [v for v in graph.vertices() if v not in support_set]
    order = graph.order()
    results = []
    for subset in Subsets(VX, j):
        test_set = support_set.union(subset)
        outcome = gzerosgame(graph, test_set, bannedset)
        if len(outcome) == order:
            results.append(frozenset(test_set))
    return results

def find_gzfs(graph, support=None, bannedset=None, upper_bound=None, lower_bound=None):
    """
    For a given graph with support and banned set, return the an optimal generalized zero forcing set. If upper_bound is less than the generalized zero forcing number then return ['wrong']. If lower_boun[...] 

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

def find_all_gzfs(graph, support=None, bannedset=None, upper_bound=None, lower_bound=None):
    """
    Return *all* optimal generalized zero forcing sets (every ZFS of minimum size).
    The minimum size is determined first using find_gZ, then gZ_leq_all enumerates
    every ZFS of that size.

    Input:
        graph: a simple graph
        support: a list of vertices of graph
        bannedset: a list of tuples representing banned edges of graph
        upper_bound: an integer, an upper bound for gZ (may shorten computation)
        lower_bound: an integer, a lower bound for gZ (may shorten computation)

    Output:
    A list of frozensets, each a minimum zero forcing set.  The list is in the
    stable iteration order produced by Subsets.  Returns an empty list only if
    no zero forcing set exists (degenerate/empty graph).

    Examples:
        sage: find_all_gzfs(graphs.PathGraph(5))
        [frozenset({0}), frozenset({4})]
        sage: find_all_gzfs(graphs.PathGraph(5), [1], [(3, 2)])
        [frozenset({0, 1, 3})]
    """
    if support is None:
        support = []
    if bannedset is None:
        bannedset = []
    min_size = find_gZ(graph, support, bannedset, upper_bound, lower_bound)
    return gZ_leq_all(graph, support, bannedset, min_size)

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
        [(('a', 0), ('b', 1), None), (('a', 1), ('b', 0), None), (('a', 1), ('b', 1), None), (('a', 1), ('b', 2), None), (('a', 2), ('b', 1), None), (('a', 2), ('b', 3), None), (('a', 3), ('b', 2), None), ([...
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
        [("a",j),("b",j)) for j in J].
        
    Examples:
        sage: J=[1,3,5];
        sage: print bridged_edges(J);
        [(('a', 1), ('b', 1)), (('a', 3), ('b', 3)), (('a', 5), ('b', 5))]
    """        
    return [(("a", j), ("b", j)) for j in J]

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

def find_Zell(g, _base=None, return_sets=False):
    """
    Return the zero forcing number of the looped graph obtained by placing
    (exactly one) loop on every vertex of g.

    Input:
        g: a simple graph, the underlying graph of the looped graph.
        _base: optional precomputed base bipartite graph (from _tilde_bipartite_base)
               for reuse across multiple calls.
        return_sets: bool (default False).
            * False  – return the zero forcing number as an integer (original behaviour).
            * True   – return a sorted list of lists.  Each inner list contains the
                       g-vertex indices (the 'a'-component) present in one minimum zero
                       forcing set of the associated tilde bipartite graph.  That is,
                       every inner list has exactly Zell(g) elements, and the outer list
                       contains every distinct such set.

    Output:
        int  when return_sets is False  – the zero forcing number Z_ell(g).
        list when return_sets is True   – a sorted list of lists of g-vertex
             indices, one list per minimum zero forcing set.

    Examples:
        sage: g = graphs.PathGraph(5)
        sage: find_Zell(g)
        1
        sage: find_Zell(g, return_sets=True)
        [[0], [4]]
        sage: g2 = Graph({0: [1], 1: [0]})
        sage: find_Zell(g2)
        1
        sage: find_Zell(g2, return_sets=True)
        [[0], [1]]
    """
    if _base is None:
        _base = _tilde_bipartite_base(g)
    V = list(g.vertices())
    h = _tilde_bipartite_with_I(_base, V)
    Yg = Y(g)
    if not return_sets:
        return find_gZ(h, Yg, []) - g.order()
    # Collect all minimum ZFS of the tilde bipartite graph (size = g.order() + Zell(g)),
    # then extract only the 'a'-components; each resulting list has size Zell(g).
    # Tilde bipartite vertices are 2-tuples: ('a', v) ∈ X(g) or ('b', v) ∈ Y(g).
    all_sets = find_all_gzfs(h, Yg, [])
    result = sorted(
        ([v[1] for v in s if v[0] == 'a'] for s in all_sets),
        key=lambda lst: tuple(sorted(lst)),
    )
    # Ensure deterministic ordering within each set
    result = [sorted(lst) for lst in result]
    return result


# ---------------------------------------------------------------------------
# Demo – run this file directly to see find_Zell in action:
#   sage diag_anal.py   (or load("diag_anal.py") inside a Sage session)
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    # Example 1: path graph on 5 vertices
    # Only the two endpoints {0} and {4} are minimum zero forcing sets.
    g_path = graphs.PathGraph(5)
    print("=== Path graph P_5 ===")
    print("find_Zell(g)              ->", find_Zell(g_path))
    print("find_Zell(g, return_sets=True) ->", find_Zell(g_path, return_sets=True))
    print()

    # Example 2: complete graph K_3
    # Every single vertex forces the rest, so all singleton sets are minima.
    g_k3 = graphs.CompleteGraph(3)
    print("=== Complete graph K_3 ===")
    print("find_Zell(g)              ->", find_Zell(g_k3))
    print("find_Zell(g, return_sets=True) ->", find_Zell(g_k3, return_sets=True))
    print()

    # Example 3: edge graph K_2 (two vertices connected by one edge)
    # Both {0} and {1} achieve the minimum.
    g_k2 = Graph({0: [1], 1: [0]})
    print("=== Edge graph K_2 ===")
    print("find_Zell(g)              ->", find_Zell(g_k2))
    print("find_Zell(g, return_sets=True) ->", find_Zell(g_k2, return_sets=True))

from itertools import permutations

def glued_graph_from_phi(G, H, phi, prefix_G="G:", prefix_H="H:"):
    VG = list(G.vertices())
    VH = list(H.vertices())

    if set(phi.keys()) != set(VG):
        raise ValueError("phi must be defined on all vertices of G.")
    if set(phi.values()) != set(VH):
        raise ValueError("phi values must be exactly the vertex set of H.")

    g = Graph()
    mapG = {u: f"{prefix_G}{u}" for u in VG}
    mapH = {v: f"{prefix_H}{v}" for v in VH}

    g.add_vertices(list(mapG.values()))
    g.add_vertices(list(mapH.values()))
    g.add_edges([(mapG[u], mapG[v]) for (u, v, _) in G.edges(labels=True)])
    g.add_edges([(mapH[u], mapH[v]) for (u, v, _) in H.edges(labels=True)])

    # glue edges u(G) -- phi(u)(H)
    g.add_edges([(mapG[u], mapH[phi[u]]) for u in VG])
    return g


def all_bijections_VG_to_VH(G, H):
    VG = list(G.vertices())
    VH = list(H.vertices())
    if len(VG) != len(VH):
        raise ValueError("Need |V(G)| = |V(H)| to have bijections.")
    for perm in permutations(VH):
        yield dict(zip(VG, perm))


def max_Zell_over_all_bijections(G, H, *, return_argmax=False):
    """
    Returns max_{bijections phi} find_Zell(glue(G,H,phi)).
    If return_argmax=True, returns (max_value, phi_that_achieves_it).
    """
    best_val = None
    best_phi = None
    # optional speed-up: reuse base bipartite for each glued graph (but glued changes each time)
    for phi in all_bijections_VG_to_VH(G, H):
        glued = glued_graph_from_phi(G, H, phi)
        val = find_Zell(glued, return_sets=False)  # integer
        if best_val is None or val > best_val:
            best_val = val
            best_phi = phi
    return (best_val, best_phi) if return_argmax else best_val

def looped_maxnull_bounds(g):
    """
    Return lower/upper bounds for looped max nullity:
      lower = max(kappa, deltaC, ncc)
      upper = Zell

    Robust to helper functions returning None.
    """
    n = g.order()
    Zell = find_Zell(g)

    kappa = g.vertex_connectivity()

    deltaC = deltaCeiling(g)
    if deltaC is None:
        deltaC = 0

    ecc = edge_clique_cover_minimum(g)
    # If no cover returned (e.g. unexpected None), fall back to weakest safe bound.
    if ecc is None:
        ncc = 0
    else:
        ncc = n - len(ecc)

    minbound = max(kappa, deltaC, ncc)
    maxbound = Zell

    return [minbound] if minbound == maxbound else [minbound, maxbound]

def edge_clique_cover_minimum(self, bound=None):
    """
    Returns an minimum edge clique cover for the graph if the
    number of covering cliques is at most ``bound``; otherwise,
    returns ``None``.

    An edge clique cover is a set of cliques which contain all of
    the edges of the graph.
    
    .. note::
        This function assumes self is connected.

    :param bound: the maximum number of cliques to consider in an
       edge clique cover

    :return: If a minimum edge clique cover is found that has at
        most ``bound`` cliques, the edge clique cover is returned
        as a list of lists, each sublist being the vertices of a
        clique. If a clique cover from this function requires more
        than ``bound`` cliques, ``None`` is returned.

    EXAMPLES::

        sage: graphs.PathGraph(3).edge_clique_cover_minimum()
        [[0, 1], [1, 2]]
        sage: graphs.CompleteGraph(5).edge_clique_cover_minimum()
        [[0, 1, 2, 3, 4]]
        sage: graphs.HouseGraph().edge_clique_cover_minimum()
        [[2, 3, 4], [0, 1], [0, 2], [1, 3]]
        sage: graphs.PetersenGraph().edge_clique_cover_minimum(bound=4)
    """
    from sage.all import ceil, Combinations
    # Take care of trivial case
    if self.size() == 0:
        return []

    max_cliques=self.cliques_maximal()
    max_cliques.sort(key=len)
    largest_clique_vertices = len(max_cliques[-1])
    max_cliques = [sorted(clique) for clique in max_cliques]
    largest_clique_edges = largest_clique_vertices \
                           *(largest_clique_vertices-1)/2
    edges_of_graph=self.edges(labels=False)
    num_edges = self.size()

    mandatory_cliques=[]


    for v in self.vertices():
        # If v is contained in only one clique, then that clique must
        # be in the clique cover
        cliques_containing_v = [c for c in max_cliques if v in c]
        if len(cliques_containing_v)==1 \
                and (cliques_containing_v[0] not in mandatory_cliques):
            mandatory_cliques.append(cliques_containing_v[0])
    for e in self.edges():
        # If e is contained in only one clique, then that clique must
        # be in the clique cover
        cliques_containing_e = [c for c in max_cliques 
                                if e[0] in c and e[1] in c]
        if len(cliques_containing_e)==1 \
                and (cliques_containing_e[0] not in mandatory_cliques):
            mandatory_cliques.append(cliques_containing_e[0])

    # Check to see if mandatory_cliques contains a clique cover
    edges_in_set_of_cliques = set([])
    for clique in mandatory_cliques:
        edges_in_clique = [(clique[i], clique[j]) 
                           for i in range(len(clique)) 
                           for j in range(i+1,len(clique))]
        edges_in_set_of_cliques.update(set(edges_in_clique))
    if len(edges_in_set_of_cliques) == num_edges:
        if bound is None or len(mandatory_cliques) <= bound:
            return mandatory_cliques
        else:
            # There are too many cliques.  Return None to be
            # consistent with the documentation, even though we
            # actually know the clique cover number (and it is greater
            # than bound).
            return None

    max_cliques = [c for c in max_cliques if c not in mandatory_cliques]
    if bound==None:
        stopping_point=len(max_cliques)
    else:
        stopping_point=min(len(max_cliques),bound-len(mandatory_cliques))

    starting_point = max(1,ceil(float(num_edges) / largest_clique_edges) \
                             - len(mandatory_cliques))
    for i in range(starting_point,stopping_point+1):
        for set_of_cliques in Combinations(max_cliques,i):
            edges_in_set_of_cliques = set([])
            for clique in set_of_cliques+mandatory_cliques:
                edges_in_clique = [(clique[i], clique[j]) 
                                   for i in range(len(clique)) 
                                   for j in range(i+1,len(clique))]
                edges_in_set_of_cliques.update(set(edges_in_clique))
            if len(edges_in_set_of_cliques) == num_edges:
                return set_of_cliques+mandatory_cliques
    return None

def contract_edge(gph,e):
    if gph.has_edge(e)==False:
        raise ValueError;
    ngh1=gph.neighbors(e[0]);
    ngh2=gph.neighbors(e[1]);
    h=gph.copy();
    h.delete_vertices([e[0],e[1]]);
    try_again=True;
    i=0;
    while try_again:
        i+=1;
        if h.has_vertex(i)==False:
            h.add_vertex(i);
            try_again=False;
    h.add_edges([(i,j) for j in set(ngh1).union(set(ngh2)).difference(set([e[0],e[1]]))]);
    return h;
    

def deltaCeiling(gph):
    gph_order=gph.order();
    delta=min(gph.degree());
    if delta==gph_order-1:
        return delta;
    name=gph.canonical_label().graph6_string();
    try:
        return dC_dic[gph_order][name];
    except: 
        pass         
    if gph.is_connected()==False:
        max_dC=delta;
        for com in gph.connected_components_subgraphs():
            max_dC=max(max_dC,deltaCeiling(com));
        return max_dC;

# ===========================================================================
# Maximum-nullity witness construction (Groebner + robust fallback)
# Drop-in block for loopedZ.py
# ===========================================================================

from itertools import combinations
from sage.all import QQ, Matrix, PolynomialRing


def _graph_pattern_symbolic_matrix_poly(G, field=QQ):
    """
    Build symbolic symmetric matrix A over polynomial ring with graph sparsity:
      A[i,j] = 0 for non-edges (i!=j)
      A[i,j] = e_{i,j} for edges (i<j, symmetric)
      A[i,i] = d_i
    Returns:
      H, R, A, diag_vars, edge_vars, vars_all
    where H is relabeled to vertices 0..n-1.
    """
    H = G.copy()
    H.relabel()
    n = H.order()

    names = [f"d_{i}" for i in range(n)]
    edge_keys = []
    for i in range(n):
        for j in range(i + 1, n):
            if H.has_edge(i, j):
                edge_keys.append((i, j))
                names.append(f"e_{i}_{j}")

    R = PolynomialRing(field, names, order="degrevlex")
    gens = list(R.gens())

    diag_vars = gens[:n]
    edge_vars = {}
    p = n
    for key in edge_keys:
        edge_vars[key] = gens[p]
        p += 1

    A = Matrix(R, n, n, 0)
    for i in range(n):
        A[i, i] = diag_vars[i]
    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) in edge_vars:
                x = edge_vars[(i, j)]
                A[i, j] = x
                A[j, i] = x
            else:
                A[i, j] = 0
                A[j, i] = 0

    return H, R, A, diag_vars, edge_vars, gens


def _all_k_minors(A, k):
    """
    All kxk minors of square matrix A.
    """
    n = A.nrows()
    if k < 0 or k > n:
        return []
    if k == 0:
        return [A.base_ring().one()]
    out = []
    idx = range(n)
    for I in combinations(idx, k):
        for J in combinations(idx, k):
            out.append(A.matrix_from_rows_and_columns(I, J).det())
    return out


def _extract_point_from_lex_groebner(Gb, R, vars_order):
    """
    Conservative extractor for lex Groebner basis:
    repeatedly solve linear univariate equations after substitution.
    Returns dict var->value or None.
    """
    sol = {}
    basis = [R(f) for f in Gb if f != 0]

    changed = True
    while changed:
        changed = False
        for f in basis:
            g = f.subs(sol) if sol else f
            if g == 0:
                continue

            rem = [v for v in vars_order if v not in sol and g.degree(v) > 0]
            if len(rem) != 1:
                continue
            v = rem[0]
            if g.degree(v) != 1:
                continue

            a = g.coefficient({v: 1})
            b = g.subs({v: 0})
            if a == 0:
                continue

            sol[v] = -b / a
            changed = True

    if len(sol) != len(vars_order):
        return None

    for f in basis:
        if f.subs(sol) != 0:
            return None
    return sol


def _repair_nonzero_diagonal_preserve_rank(M, max_tries=40, step_values=None):
    """
    Try to make every diagonal entry nonzero by diagonal perturbations while
    preserving rank(M). Returns repaired matrix or None.
    """
    if step_values is None:
        step_values = [1, -1, 2, -2, 3, -3, QQ(1, 2), QQ(-1, 2), QQ(3, 2), QQ(-3, 2)]

    A = Matrix(QQ, M)
    target_rank = A.rank()
    n = A.nrows()

    if all(A[i, i] != 0 for i in range(n)):
        return A

    for _ in range(max_tries):
        changed = False
        for i in range(n):
            if A[i, i] != 0:
                continue
            fixed = False
            for t in step_values:
                B = Matrix(QQ, A)
                B[i, i] += QQ(t)
                if B[i, i] != 0 and B.rank() == target_rank:
                    A = B
                    changed = True
                    fixed = True
                    break
            if not fixed:
                pass

        if all(A[i, i] != 0 for i in range(n)):
            return A
        if not changed:
            break

    return None


def max_nullity_witness_matrix_groebner(
    G,
    target_nullity=None,
    require_nz_diag=False,
    field=QQ,
    random_fallback_trials=3000,
    random_range=(-5, 5),
    return_diagnostics=False,
):
    """
    Construct a symmetric matrix matching G's pattern and (attempt to) maximize nullity.

    Constraints:
      - offdiag non-edge => 0
      - offdiag edge     => nonzero
      - diagonal free; if require_nz_diag=True then diagonal must be nonzero

    Returns dict:
      {
        'matrix': Matrix or None,
        'rank': int or None,
        'nullity': int or None,
        'assignment': dict or None,
        'method': 'groebner',
        'diagnostics': ... (optional)
      }
    """
    H, R0, A0, diag_vars0, edge_vars0, vars0 = _graph_pattern_symbolic_matrix_poly(G, field=field)
    n = H.order()

    diagnostics = {"attempts": []}
    best = {"matrix": None, "rank": None, "nullity": -1, "assignment": None, "method": "groebner"}

    # Search low rank first => high nullity first
    for r in range(0, n + 1):
        k = r + 1
        if k > n:
            continue

        minors0 = _all_k_minors(A0, k)

        # Lex ring for elimination/extraction
        R = PolynomialRing(field, [str(v) for v in vars0], order="lex")
        gens = list(R.gens())
        vmap = {str(vars0[i]): gens[i] for i in range(len(vars0))}
        sub_map = {vars0[i]: vmap[str(vars0[i])] for i in range(len(vars0))}

        A = Matrix(R, n, n, lambda i, j: R(A0[i, j].subs(sub_map)))
        eqs = [R(m.subs(sub_map)) for m in minors0]

        I = R.ideal(eqs)
        att = {"rank_bound": r, "num_eqs": len(eqs)}

        try:
            Glex = I.groebner_basis()
            att["groebner_ok"] = True
            att["groebner_len"] = len(Glex)
        except Exception as e:
            att["groebner_ok"] = False
            att["groebner_error"] = str(e)
            diagnostics["attempts"].append(att)
            continue

        if R.one() in Glex:
            att["inconsistent"] = True
            diagnostics["attempts"].append(att)
            continue

        sol = _extract_point_from_lex_groebner(Glex, R, gens)

        # Robust randomized fallback on minors equations
        if sol is None and random_fallback_trials > 0:
            import random
            lo, hi = random_range
            nz_choices = [x for x in range(lo, hi + 1) if x != 0] or [-1, 1]

            edge_names = {f"e_{i}_{j}" for (i, j) in edge_vars0.keys()}
            diag_names = {f"d_{i}" for i in range(n)}

            for _ in range(random_fallback_trials):
                cand = {}
                for v in gens:
                    sv = str(v)
                    if sv in edge_names or (require_nz_diag and sv in diag_names):
                        cand[v] = random.choice(nz_choices)
                    else:
                        cand[v] = random.randint(lo, hi)

                ok = True
                for f in eqs:
                    if f.subs(cand) != 0:
                        ok = False
                        break
                if ok:
                    sol = cand
                    att["used_random_fallback"] = True
                    break

        if sol is None:
            att["solution_found"] = False
            diagnostics["attempts"].append(att)
            continue

        M = Matrix(QQ, n, n, lambda i, j: QQ(A[i, j].subs(sol)))

        # Repair diagonal if requested
        if require_nz_diag and any(M[i, i] == 0 for i in range(n)):
            repaired = _repair_nonzero_diagonal_preserve_rank(M)
            if repaired is not None:
                M = repaired

        # Validate pattern
        valid = True
        for i in range(n):
            if require_nz_diag and M[i, i] == 0:
                valid = False
                break
            for j in range(i + 1, n):
                if H.has_edge(i, j):
                    if M[i, j] == 0:
                        valid = False
                        break
                else:
                    if M[i, j] != 0:
                        valid = False
                        break
            if not valid:
                break

        rr = int(M.rank())
        nn = n - rr

        att["solution_found"] = True
        att["valid_pattern"] = valid
        att["rank"] = rr
        att["nullity"] = nn
        diagnostics["attempts"].append(att)

        if valid and nn > best["nullity"]:
            best = {"matrix": M, "rank": rr, "nullity": nn, "assignment": sol, "method": "groebner"}

        # early stop if target reached
        if valid and target_nullity is not None and nn >= target_nullity:
            out = best
            if return_diagnostics:
                out = dict(out)
                out["diagnostics"] = diagnostics
            return out

        # if this rank bound is achieved by a valid point, that's usually near-optimal in this search order
        if valid and target_nullity is None:
            out = best
            if return_diagnostics:
                out = dict(out)
                out["diagnostics"] = diagnostics
            return out

    if best["matrix"] is None:
        out = {"matrix": None, "rank": None, "nullity": None, "assignment": None, "method": "groebner"}
    else:
        out = best

    if return_diagnostics:
        out = dict(out)
        out["diagnostics"] = diagnostics
    return out    
    max_dC=delta;
    for e in gph.edges():
        max_dC=max(max_dC,deltaCeiling(contract_edge(gph,e)));
    return max_dC;    

def max_nullity_witness_close_to_Zell(
    G,
    require_nz_diag=False,
    zell_target=None,
    max_drop=6,                 # try target, target-1, ..., target-max_drop
    attempts_per_level=8,       # repeated solver restarts per target level
    initial_trials=2000,        # fallback trials inside each attempt
    trial_growth=2,
    random_range=(-6, 6),
    return_diagnostics=False,
):
    """
    Try hard to find a witness matrix with nullity as close as possible to find_Zell(G).
    Returns best witness found, preferring exact target match.
    """
    n = G.order()
    target = find_Zell(G) if zell_target is None else int(zell_target)

    best = {
        "matrix": None,
        "rank": None,
        "nullity": -1,
        "assignment": None,
        "method": "groebner-close-wrapper",
    }
    diag = {"target": target, "levels": []}

    # desired nullities: target downwards
    for drop in range(0, max_drop + 1):
        desired_nullity = target - drop
        if desired_nullity < 0:
            break

        level_info = {"desired_nullity": desired_nullity, "attempts": []}
        trials = int(initial_trials)

        for a in range(attempts_per_level):
            res = max_nullity_witness_matrix_groebner(
                G,
                target_nullity=desired_nullity,  # early exit if this level is achieved
                require_nz_diag=require_nz_diag,
                random_fallback_trials=trials,
                random_range=random_range,
                return_diagnostics=False,
            )

            got = res.get("nullity", None)
            level_info["attempts"].append({
                "attempt": a + 1,
                "trials": trials,
                "got_nullity": got,
                "got_matrix": res.get("matrix", None) is not None,
            })

            if res.get("matrix", None) is not None and got is not None:
                got = int(got)
                if got > best["nullity"]:
                    best = {
                        "matrix": res["matrix"],
                        "rank": int(res["rank"]),
                        "nullity": got,
                        "assignment": res.get("assignment"),
                        "method": "groebner-close-wrapper",
                    }

                # exact hit on Zell: done
                if got == target:
                    out = dict(best)
                    out["target_nullity"] = target
                    out["gap_to_target"] = target - best["nullity"]
                    if return_diagnostics:
                        diag["levels"].append(level_info)
                        out["diagnostics"] = diag
                    return out

            trials *= int(trial_growth)

        diag["levels"].append(level_info)

    # final best
    if best["matrix"] is None:
        out = {
            "matrix": None,
            "rank": None,
            "nullity": None,
            "assignment": None,
            "method": "groebner-close-wrapper",
            "target_nullity": target,
            "gap_to_target": None,
        }
    else:
        out = dict(best)
        out["target_nullity"] = target
        out["gap_to_target"] = target - best["nullity"]

    if return_diagnostics:
        out["diagnostics"] = diag
    return out

def max_nullity_witness_factorized(
    G,
    target_nullity=None,
    require_nz_diag=False,
    max_drop=8,
    attempts_per_rank=30,
    random_range=(-5, 5),
    seed=None,
    return_diagnostics=False,
):
    r"""
    Find a symmetric matrix A with graph pattern via factorization A = U*D*U^T.

    Pattern:
      - A[i,j] = 0 for i!=j and non-edge
      - A[i,j] != 0 for i!=j and edge
      - if require_nz_diag: A[i,i] != 0

    Optimization objective:
      get nullity(A) as close as possible to target_nullity (default find_Zell(G))
      by trying low ranks r = n-target, n-(target-1), ... (increasing as needed).

    Returns dict:
      matrix, rank, nullity, target_nullity, gap_to_target, method, diagnostics(optional)
    """
    import random
    from sage.all import QQ, Matrix

    if seed is not None:
        random.seed(seed)

    H = G.copy()
    H.relabel()
    n = H.order()

    target = find_Zell(H) if target_nullity is None else int(target_nullity)
    lo, hi = random_range
    vals = [x for x in range(lo, hi + 1) if x != 0]
    if not vals:
        vals = [-1, 1]

    # Build edge/non-edge index lists
    edge_pairs = []
    nonedge_pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            if H.has_edge(i, j):
                edge_pairs.append((i, j))
            else:
                nonedge_pairs.append((i, j))

    best = {
        "matrix": None,
        "rank": None,
        "nullity": -1,
        "target_nullity": target,
        "gap_to_target": None,
        "method": "factorized",
    }
    diag = {"trials": []}

    # desired nullity target, then relax downward
    for drop in range(max_drop + 1):
        desired_nullity = target - drop
        if desired_nullity < 0:
            break
        r = n - desired_nullity
        if r < 0 or r > n:
            continue

        for t in range(attempts_per_rank):
            # Random U (n x r), D diagonal invertible-ish (r x r)
            if r == 0:
                # only zero matrix, useful only for edgeless graph
                A = Matrix(QQ, n, n, 0)
            else:
                U = Matrix(QQ, n, r, lambda i, j: random.choice(vals))
                Ddiag = [QQ(random.choice(vals)) for _ in range(r)]
                D = Matrix(QQ, r, r, 0)
                for k in range(r):
                    D[k, k] = Ddiag[k]
                A = U * D * U.transpose()

            # Enforce zero on non-edges by rejection
            ok = True
            for (i, j) in nonedge_pairs:
                if A[i, j] != 0:
                    ok = False
                    break
            if not ok:
                if return_diagnostics:
                    diag["trials"].append({"rank_target": r, "attempt": t + 1, "status": "reject_nonedge"})
                continue

            # Enforce nonzero on edges
            for (i, j) in edge_pairs:
                if A[i, j] == 0:
                    ok = False
                    break
            if not ok:
                if return_diagnostics:
                    diag["trials"].append({"rank_target": r, "attempt": t + 1, "status": "reject_edge_zero"})
                continue

            # Diagonal nonzero condition if requested
            if require_nz_diag:
                if any(A[i, i] == 0 for i in range(n)):
                    # try your existing repair if available
                    try:
                        repaired = _repair_nonzero_diagonal_preserve_rank(A)
                    except Exception:
                        repaired = None
                    if repaired is None or any(repaired[i, i] == 0 for i in range(n)):
                        if return_diagnostics:
                            diag["trials"].append({"rank_target": r, "attempt": t + 1, "status": "reject_diag_zero"})
                        continue
                    A = repaired

            rr = int(A.rank())
            nn = n - rr

            if return_diagnostics:
                diag["trials"].append({
                    "rank_target": r, "attempt": t + 1, "status": "accept",
                    "rank": rr, "nullity": nn
                })

            if nn > best["nullity"]:
                best.update({
                    "matrix": A,
                    "rank": rr,
                    "nullity": nn,
                    "gap_to_target": target - nn,
                })

            if nn >= target:
                if return_diagnostics:
                    best["diagnostics"] = diag
                return best

    if return_diagnostics:
        best["diagnostics"] = diag
    return best
