"""Discoid graphs: definitions and exact test.
0-dominated v: exists w in N(v) with N[v] subset N[w]   (<=> G[N(v)] is a cone)
0-dismantlable: reduces to K_1 by removing 0-dominated vertices
1-dominated v: N(v) nonempty and G[N(v)] is 0-dismantlable
discoid (= 1-dismantlable): reduces to K_1 by removing 1-dominated vertices (some order)
"""
import itertools, sys, random
import networkx as nx

def to_adj(G):
    nodes = list(G.nodes()); idx = {v: i for i, v in enumerate(nodes)}
    adj = [0] * len(nodes)
    for u, v in G.edges():
        adj[idx[u]] |= 1 << idx[v]; adj[idx[v]] |= 1 << idx[u]
    return adj, nodes

def popcount(x): return bin(x).count('1')

def zero_dominated(adj, S, v):
    """v in S is 0-dominated in G[S]"""
    Nv = (adj[v] | (1 << v)) & S
    m = adj[v] & S
    while m:
        b = m & -m; w = b.bit_length() - 1; m ^= b
        if Nv & ~((adj[w] | b) & S) == 0: return True
    return False

def zero_dismantlable(adj, S):
    """G[S] reduces to K_1 by removing 0-dominated vertices (order irrelevant: the core is unique)"""
    if S == 0: return False
    changed = True
    while changed and popcount(S) > 1:
        changed = False
        m = S
        while m:
            b = m & -m; v = b.bit_length() - 1; m ^= b
            if zero_dominated(adj, S, v):
                S &= ~b; changed = True; break
    return popcount(S) == 1

def one_dominated(adj, S, v):
    N = adj[v] & S
    return N != 0 and zero_dismantlable(adj, N)

class Discoid:
    def __init__(self, adj, cap=400000):
        self.adj = adj; self.memo = {}; self.cap = cap; self.aborted = False
    def test(self, S, path=None):
        """returns True/False, or None if aborted (state cap)."""
        if popcount(S) == 1: return True
        if S in self.memo: return self.memo[S]
        if len(self.memo) > self.cap: self.aborted = True; return None
        # candidates: 1-dominated vertices; try 0-dominated first (cheap, always safe? not nec.)
        cands = []
        m = S
        while m:
            b = m & -m; v = b.bit_length() - 1; m ^= b
            if one_dominated(self.adj, S, v): cands.append((0 if zero_dominated(self.adj, S, v) else 1, v))
        cands.sort()
        res = False
        for _, v in cands:
            r = self.test(S & ~(1 << v))
            if r is None: res = None; break
            if r:
                res = True
                if path is not None: path.append(v)
                break
        self.memo[S] = res
        return res

def is_discoid(G, cap=400000, with_order=False):
    adj, nodes = to_adj(G)
    D = Discoid(adj, cap)
    S = (1 << len(adj)) - 1
    r = D.test(S)
    if not with_order: return r
    if not r: return r, None
    # reconstruct an order greedily using memo
    order = []
    while popcount(S) > 1:
        for v in range(len(adj)):
            if (S >> v) & 1 and one_dominated(adj, S, v) and D.memo.get(S & ~(1 << v), popcount(S & ~(1 << v)) == 1):
                order.append(nodes[v]); S &= ~(1 << v); break
    return r, order

def clique_graph(G):
    C = [frozenset(c) for c in nx.find_cliques(G)]
    K = nx.Graph(); K.add_nodes_from(range(len(C)))
    for i, j in itertools.combinations(range(len(C)), 2):
        if C[i] & C[j]: K.add_edge(i, j)
    return K, C
