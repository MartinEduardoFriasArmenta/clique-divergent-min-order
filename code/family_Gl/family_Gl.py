"""Verification of the family G_l (l >= 1) of Section 6 ("A family on a flag 2-sphere with three apices").

T_l: flag 2-sphere formed by 4-cycles C^0,...,C^l (vertex (k,i) = vertex i of C^k), consecutive cycles joined by
antiprism bands ((k,i)~(k+1,i), (k,i)~(k+1,i+1)), and caps N (joined to C^0) and S (joined to C^l).
G_l = T_l + {a,b,c}: a,b,c pairwise nonadjacent; c ~ all of T_l; a ~ T_l - {x} with x=(0,0); b ~ T_l - {N}.
|V(G_l)| = 4l+9, |E(G_l)| = 24l+28, |K(G_l)| = 24l+15.

Checks:
 (1) the order N, C^0, ..., C^l, a, b, S deletes 1-dominated vertices and leaves c            (l = 1..L)
 (2) T_l is a flag 2-sphere; G_l has no 0-dominated vertex, is not a suspension, and has no clique separator (l = 1..10)
 (3) the explicit map r_l (identity on D_l = C^{l-1} u C^l u {S} and on b, c; a, N -> c; (l-2,i) -> (l,i+1);
     (k,i) -> S for k <= l-3) is a retraction of G_l onto G_l[D_l + {b,c}] = Susp(D_l), which is isomorphic to G_LPV
     (graph6 JrqipkN~~~?)                                                                          (l = 1..L)
 (4) independent check: K(G_l) retracts onto an induced O_4:
       l = 1..8 : explicit retraction found by constraint satisfaction and checked edge by edge;
       l >= 9   : uniform retraction rho_l = fixed table on the "window" cliques (tube vertices in layers <= 7,
                  N counting as layer -1), far cliques sent to ALPHA (apex a or b) or GAMMA (apex c).
                  For l = 9..L we check that rho_l is a retraction and that the set of edges of K(G_l) having an
                  endpoint in the window, together with their images, is the same as for l = 9 (cliques with no
                  vertex in layers <= 8 are described only by their apex). Since edges between far cliques are sent
                  to {ALPHA, GAMMA}, which are adjacent, this proves the statement for every l >= 9.
 Certificates are written to family_Gl_certificates.json.
Usage: python3 family_Gl.py [L]      (default L = 16; runs in a few seconds)
"""
import itertools, sys, json
import networkx as nx
sys.path.insert(0, '.')
from discoid import to_adj, one_dominated, zero_dominated, popcount
from retract_oct import retract_octahedron, retract_octahedron_fixed

APICES = ('a', 'b', 'c')
X = (0, 0)

def tube(l, m=4):
    T = nx.Graph()
    for k in range(l + 1):
        for i in range(m): T.add_edge((k, i), (k, (i + 1) % m))
    for k in range(l):
        for i in range(m): T.add_edge((k, i), (k + 1, i)); T.add_edge((k, i), (k + 1, (i + 1) % m))
    for i in range(m): T.add_edge('N', (0, i)); T.add_edge('S', (l, i))
    return T

def G_(l):
    T = tube(l); G = nx.Graph(); G.add_edges_from(T.edges())
    for v in T:
        G.add_edge('c', v)
        if v != X: G.add_edge('a', v)
        if v != 'N': G.add_edge('b', v)
    return G

def is_flag_sphere(T):
    n, e = T.number_of_nodes(), T.number_of_edges()
    cl = list(nx.enumerate_all_cliques(T))
    return nx.check_planarity(T)[0] and e == 3 * n - 6 and sum(len(c) == 3 for c in cl) == 2 * n - 4 and not any(len(c) == 4 for c in cl)

def has_clique_separator(G):
    for C in nx.enumerate_all_cliques(G):
        H = G.copy(); H.remove_nodes_from(C)
        if H.number_of_nodes() and not nx.is_connected(H): return True
    return False

def layer(v): return -1 if v == 'N' else (10 ** 6 if v == 'S' else v[0])

def cliques(G):
    out = []
    for c in nx.find_cliques(G):
        ap = [v for v in c if v in APICES]
        assert len(ap) == 1 and len(c) == 4, c            # every clique = one apex + a face of T_l
        out.append((ap[0], frozenset(v for v in c if v not in APICES)))
    return out

def Kgraph(C):
    K = nx.Graph(); K.add_nodes_from(C)
    for P, Q in itertools.combinations(C, 2):
        if P[0] == Q[0] or P[1] & Q[1]: K.add_edge(P, Q)
    return K

def deletion_order(l): return ['N'] + [(k, i) for k in range(l + 1) for i in range(4)] + ['a', 'b', 'S']

def check_order(l):
    G = G_(l); adj, nodes = to_adj(G); idx = {v: i for i, v in enumerate(nodes)}; S = (1 << len(nodes)) - 1
    for v in deletion_order(l):
        if not one_dominated(adj, S, idx[v]): return False
        S &= ~(1 << idx[v])
    return popcount(S) == 1 and nodes[S.bit_length() - 1] == 'c'

def basic(l):
    G = G_(l); n = G.number_of_nodes(); adj, nodes = to_adj(G); full = (1 << n) - 1
    dom = sum(1 for v in range(n) if zero_dominated(adj, full, v))
    susp = any(G.degree(u) == n - 2 == G.degree(w) and not G.has_edge(u, w) for u, w in itertools.combinations(G, 2))
    return n, G.number_of_edges(), is_flag_sphere(tube(l)), dom, susp, has_clique_separator(G)

def is_retraction(K, S, r):
    anti = {S[i]: S[i ^ 1] for i in range(len(S))}
    induced = all((not K.has_edge(S[i], S[j])) == (j == i ^ 1) for i in range(len(S)) for j in range(len(S)) if i != j)
    return induced and all(r[s] == s for s in S) and all(r[u] != anti[r[v]] for u, v in K.edges())

def octahedra(K, k):
    V = sorted(K, key=str); ne = [(u, v) for u, v in itertools.combinations(V, 2) if not K.has_edge(u, v)]
    def rec(s, ch):
        if len(ch) == 2 * k: yield list(ch); return
        for i in range(s, len(ne)):
            u, v = ne[i]
            if u in ch or v in ch: continue
            if all(K.has_edge(u, y) and K.has_edge(v, y) for y in ch): yield from rec(i + 1, ch + [u, v])
    yield from rec(0, [])

def csp_certificate(l):
    G = G_(l); C = cliques(G); K = Kgraph(C)
    for Sx in octahedra(K, 4):
        r = retract_octahedron(K, Sx, limit=300000)
        if r and r != 'limite' and is_retraction(K, Sx, r): return Sx, r
    return None, None

def r_explicit(l):
    r = {}
    for v in G_(l):
        if v in ('b', 'c', 'S') or (isinstance(v, tuple) and v[0] >= l - 1): r[v] = v
        elif v in ('a', 'N'): r[v] = 'c'
        elif v[0] == l - 2: r[v] = (l, (v[1] + 1) % 4)
        else: r[v] = 'S'
    return r

def check_explicit(l):
    G = G_(l); r = r_explicit(l); H = G.subgraph(set(r.values()))
    morph = all(r[u] == r[v] or H.has_edge(r[u], r[v]) for u, v in G.edges())
    ident = all(r[v] == v for v in H)
    return morph and ident and nx.is_isomorphic(H, nx.from_graph6_bytes(b'JrqipkN~~~?'))

def window(Q): return all(layer(v) <= 7 for v in Q[1])
def enc(Q): return [Q[0], sorted(map(str, Q[1]))]

if __name__ == '__main__':
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    cert = {'small': {}, 'uniform': None}
    print('(1) deletion order N, C^0, ..., C^l, a, b, S of 1-dominated vertices leaving c:',
          all(check_order(l) for l in range(1, L + 1)), '(l = 1..%d)' % L, flush=True)
    for l in range(1, 11):
        n, e, flag, dom, susp, sep = basic(l)
        print('(2) l=%2d n=%2d m=%3d flag sphere=%s 0-dominated=%d suspension=%s clique separator=%s' % (l, n, e, flag, dom, susp, sep), flush=True)
        assert n == 4 * l + 9 and e == 24 * l + 28 and flag and dom == 0 and not susp and not sep
    print('(3) explicit retraction G_l -> Susp(D_l) = G_LPV:', all(check_explicit(l) for l in range(1, L + 1)), '(l = 1..%d)' % L, flush=True)
    for l in range(1, 9):
        Sx, r = csp_certificate(l)
        print('(4) l=%d: |K(G_l)|=%d, retraction K(G_l) -> O_4 found and checked: %s' % (l, 24 * l + 15, Sx is not None), flush=True)
        cert['small'][l] = {'O4': [enc(Q) for Q in Sx], 'r': [[enc(Q), enc(P)] for Q, P in r.items()]}
    # uniform retraction for l >= 9
    G9 = G_(9); C9 = cliques(G9); K9 = Kgraph(C9)
    q = lambda ap, *vs: (ap, frozenset(vs))
    Sx = [q('a', (4, 0), (4, 1), (5, 1)), q('c', (3, 1), (3, 2), (4, 2)),
          q('a', (3, 0), (3, 1), (4, 1)), q('c', (3, 3), (4, 0), (4, 3)),
          q('a', (4, 2), (4, 3), (5, 3)), q('c', (3, 0), (3, 3), (4, 0)),
          q('a', (3, 2), (3, 3), (4, 3)), q('c', (3, 1), (4, 1), (4, 2))]      # antipodal pairs (0,1),(2,3),(4,5),(6,7)
    ALPHA, GAMMA = Sx[0], Sx[5]
    fixed = {Q: (ALPHA if Q[0] in 'ab' else GAMMA) for Q in C9 if not window(Q) and Q not in Sx}
    r9 = retract_octahedron_fixed(K9, Sx, fixed)
    assert r9 and r9 != 'limite' and is_retraction(K9, Sx, r9)
    table = {Q: r9[Q] for Q in C9 if window(Q)}
    print('    window cliques: %d; ALPHA = %s, GAMMA = %s, adjacent: %s' % (len(table), enc(ALPHA), enc(GAMMA), K9.has_edge(ALPHA, GAMMA)), flush=True)
    ref = None
    for l in range(9, L + 1):
        G = G_(l); C = cliques(G); K = Kgraph(C)
        rho = {Q: (table[Q] if window(Q) else (ALPHA if Q[0] in 'ab' else GAMMA)) for Q in C}
        ok = is_retraction(K, Sx, rho)
        desc = lambda Q: Q if any(layer(v) <= 8 for v in Q[1]) else ('far', Q[0])
        loc = set()
        for P, Q in K.edges():
            if window(P): loc.add((P, rho[P], desc(Q), rho[Q]))
            if window(Q): loc.add((Q, rho[Q], desc(P), rho[P]))
        loc = frozenset(loc); same = ref is None or loc == ref
        if ref is None: ref = loc
        print('(4) l=%2d |K(G_l)|=%d  rho_l is a retraction onto O_4: %s ; local edges identical to l=9: %s' % (l, len(C), ok, same), flush=True)
        assert ok and same and len(C) == 24 * l + 15
    cert['uniform'] = {'O4': [enc(Q) for Q in Sx], 'ALPHA': enc(ALPHA), 'GAMMA': enc(GAMMA),
                       'window_table': [[enc(Q), enc(P)] for Q, P in table.items()]}
    json.dump(cert, open('family_Gl_certificates.json', 'w'))
    print('all checks passed; certificates written to family_Gl_certificates.json')
