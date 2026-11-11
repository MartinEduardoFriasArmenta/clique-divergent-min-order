import networkx as nx, itertools, sys
import numpy as np
def betti_z2(G, maxdim=4):
    # simplices = complete subgraphs up to dim maxdim
    simp={}
    for k in range(1,maxdim+2):
        simp[k-1]=[frozenset(c) for c in nx.enumerate_all_cliques(G) if len(c)==k]
    def rank_z2(M):
        M=M.copy()%2; r=0; rows,cols=M.shape
        for c in range(cols):
            piv=None
            for i in range(r,rows):
                if M[i,c]: piv=i; break
            if piv is None: continue
            M[[r,piv]]=M[[piv,r]]
            for i in range(rows):
                if i!=r and M[i,c]: M[i]^=M[r]
            r+=1
            if r==rows: break
        return r
    ranks={}
    idx={d:{s:i for i,s in enumerate(simp[d])} for d in simp}
    for d in range(1,maxdim+1):
        if not simp[d] or not simp[d-1]: ranks[d]=0; continue
        M=np.zeros((len(simp[d-1]),len(simp[d])),dtype=np.uint8)
        for j,s in enumerate(simp[d]):
            for v in s:
                M[idx[d-1][s-{v}],j]=1
        ranks[d]=rank_z2(M)
    b=[]
    for d in range(0,maxdim+1):
        n=len(simp[d]); rk=ranks.get(d,0); rk1=ranks.get(d+1,0)
        b.append(n-rk-rk1)
    b[0]-=1  # reduced
    return b
if __name__=="__main__":
    # usage: python3 homol.py <graph6>  -> reduced Betti numbers over Z_2 of the clique complex
    G=nx.from_graph6_bytes(sys.argv[1].encode())
    print(betti_z2(G))
