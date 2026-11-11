import networkx as nx, itertools, sys
from discoid import clique_graph
def octahedra(G, k):
    """induced K_{2,...,2} with k parts: sets of 2k vertices whose induced graph is complete minus a perfect matching"""
    res=[]
    for S in itertools.combinations(sorted(G),2*k):
        H=G.subgraph(S)
        if H.number_of_edges()!=2*k*(k-1) or any(d!=2*k-2 for _,d in H.degree()): continue
        res.append(S)
    return res
def retraction(G,S):
    """find r: V(G)->S, identity on S, with uv edge => r(u)=r(v) or r(u)r(v) edge. backtracking"""
    S=list(S); others=[v for v in G if v not in S]; r={v:v for v in S}
    def ok(v,x):
        for u in G[v]:
            if u in r and r[u]!=x and not G.has_edge(r[u],x): return False
        return True
    def bt(i):
        if i==len(others): return True
        v=others[i]
        for x in S:
            if ok(v,x):
                r[v]=x
                if bt(i+1): return True
                del r[v]
        return False
    return dict(r) if bt(0) else None
if __name__=="__main__":
    A=nx.from_graph6_bytes(b'J~Z}eY~qzd_')
    for name,G in (("A",A),("K(A)",clique_graph(A)[0])):
        for k in (3,4):
            octs=octahedra(G,k); print(name,"octaedros inducidos O_%d:"%k,len(octs),flush=True)
            for S in octs:
                r=retraction(G,S)
                if r: print("  RETRACTO:",S,{v:x for v,x in r.items() if v!=x}); 
