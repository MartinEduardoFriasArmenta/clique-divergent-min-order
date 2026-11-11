"""Experimento: para triangulaciones de Whitney del disco B, ¿son B, K(B), ..., K^5(B) 1-desmantelables?
Uso: python3 disk_exp.py <geng> <minutos> <salida.tsv>
- Discos con 3..9 vértices: todos (geng + prueba exacta de que Delta(B) es un disco).
- Discos con 10..30 vértices: generados al azar por 'shelling' inverso (random_disk), hasta agotar el tiempo.
Prueba de 1-desmantelabilidad: borrado voraz aleatorio (certificado positivo) y, si falla, búsqueda exacta con memo.
"""
import sys, time, subprocess, random, itertools
import networkx as nx
from exp15 import parse_g6
from discoid import to_adj, zero_dominated, zero_dismantlable, one_dominated, Discoid, popcount, clique_graph
import random

def random_disk(nsteps, seed, p_ear=0.5):
    rnd=random.Random(seed)
    G=nx.Graph(); G.add_edges_from([(0,1),(1,2),(0,2)]); bd=[0,1,2]; nxt=3
    for _ in range(nsteps):
        L=len(bd)
        if rnd.random()<p_ear and L>=4:
            # try ear at a random position: bd[i-1], bd[i], bd[i+1]
            order=list(range(L)); rnd.shuffle(order); done=False
            for i in order:
                a,f,b=bd[i-1],bd[i],bd[(i+1)%L]
                if G.has_edge(a,b): continue
                if set(G[a])&set(G[b])!={f}: continue
                G.add_edge(a,b); bd.pop(i); done=True; break
            if done: continue
        i=rnd.randrange(L); a,b=bd[i],bd[(i+1)%L]
        G.add_node(nxt); G.add_edge(nxt,a); G.add_edge(nxt,b); bd.insert(i+1,nxt); nxt+=1
    return G, bd


def is_disk(G):
    """Delta(G) es una triangulación (de Whitney) del disco."""
    if G.number_of_nodes() < 3 or not nx.is_connected(G): return False
    C = [frozenset(c) for c in nx.find_cliques(G)]
    if any(len(c) > 3 for c in C): return False
    T = [c for c in C if len(c) == 3]
    tri_of_edge = {}
    for t in T:
        for e in itertools.combinations(sorted(t), 2): tri_of_edge.setdefault(e, []).append(t)
    for e in G.edges():
        k = len(tri_of_edge.get(tuple(sorted(e)), []))
        if k == 0 or k > 2: return False
    bd = nx.Graph([e for e in G.edges() if len(tri_of_edge[tuple(sorted(e))]) == 1])
    if bd.number_of_edges() == 0: return False
    for v in G:
        L = G.subgraph(G[v])
        if not nx.is_connected(L): return False
        degs = sorted(d for _, d in L.degree())
        if not (all(d == 2 for d in degs) or (degs[:2] == [1, 1] and all(d == 2 for d in degs[2:]))): return False
    if G.number_of_nodes() - G.number_of_edges() + len(T) != 1: return False
    return nx.is_connected(bd) and all(d == 2 for _, d in bd.degree())

def core0(G):
    """0-núcleo: borrado repetido de vértices 0-dominados (único salvo isomorfismo). El 0-núcleo de un disco
    es un punto o un disco cuyos vértices de la frontera tienen grado >= 4 (sin vértices 0-dominados)."""
    G = G.copy(); ch = True
    while ch:
        ch = False
        for v in list(G):
            Nv = set(G[v]) | {v}
            if any(Nv <= set(G[w]) | {w} for w in G[v]): G.remove_node(v); ch = True; break
    return G

def greedy_1dism(adj, S, rnd, tries=6):
    for t in range(tries):
        T = S
        while popcount(T) > 1:
            c0 = []; c1 = []; m = T
            while m:
                b = m & -m; v = b.bit_length() - 1; m ^= b
                if zero_dominated(adj, T, v): c0.append(v)
                elif one_dominated(adj, T, v): c1.append(v)
            pool = c0 if c0 else c1
            if not pool: break
            T &= ~(1 << rnd.choice(pool))
        if popcount(T) == 1: return True
    return None

def one_dism(G, rnd, cap=100000):
    adj, _ = to_adj(G); n = len(adj); S = (1 << n) - 1
    if zero_dismantlable(adj, S): return '0'      # 0-desmantelable (luego 1-desmantelable)
    if greedy_1dism(adj, S, rnd): return '1'      # 1-desmantelable (certificado voraz)
    D = Discoid(adj, cap=cap); r = D.test(S)
    return '1' if r else ('NO' if r is False else '?')

def analyze(B, rnd, maxit=5, cap_vertices=400):
    row = []; G = B
    for i in range(maxit + 1):
        n = G.number_of_nodes()
        if n == 1: row.append('K1'); break
        if n > cap_vertices: row.append('>%d' % cap_vertices); break
        row.append('%d:%s' % (n, one_dism(G, rnd)))
        if i < maxit:
            C = []
            for c in nx.find_cliques(G):          # enumeración con tope: si hay más de cap_vertices clanes, paramos
                C.append(frozenset(c))
                if len(C) > cap_vertices: break
            if len(C) > cap_vertices: row.append('>%d' % cap_vertices); break
            K = nx.Graph(); K.add_nodes_from(range(len(C)))
            K.add_edges_from((a, b) for a, b in itertools.combinations(range(len(C)), 2) if C[a] & C[b])
            G = K
    return row

if __name__ == '__main__':
    geng = sys.argv[1]; deadline = time.time() + 60 * float(sys.argv[2])
    resume = len(sys.argv) > 4            # python3 disk_exp.py geng minutos salida [semilla]: continúa sólo con discos aleatorios
    out = open(sys.argv[3], 'a' if resume else 'w')
    rnd = random.Random(int(sys.argv[4]) if resume else 1); count = 0
    def emit(src, B):
        global count
        row = analyze(B, rnd); count += 1
        out.write('%s\t%d\t%d\t%s\t%s\n' % (src, B.number_of_nodes(), B.number_of_edges(), nx.to_graph6_bytes(B, header=False).strip().decode(), '\t'.join(row))); out.flush()
        print(src, B.number_of_nodes(), row, flush=True)
    for n in (range(3, 10) if not resume else []):
        p = subprocess.Popen([geng, '-c', '-q', str(n)], stdout=subprocess.PIPE, text=True)
        for line in p.stdout:
            if time.time() > deadline: break
            adj = parse_g6(line); G = nx.Graph(); G.add_nodes_from(range(n))
            for i in range(n):
                for j in range(i + 1, n):
                    if (adj[i] >> j) & 1: G.add_edge(i, j)
            if G.number_of_edges() > 3 * n - 6 or not is_disk(G): continue
            emit('todos n=%d' % n, G)
        p.kill()
        if time.time() > deadline: break
    seen = set(); k = 0 if not resume else 100000 * int(sys.argv[4])
    while time.time() < deadline:
        steps = rnd.randint(5, 60); G, bd = random_disk(steps, 1000 + k); k += 1
        hard = (k % 3 != 0)
        if hard:
            G = core0(G)                      # disco sin vértices 0-dominados (o trivial)
            if G.number_of_nodes() < 8: continue
        assert is_disk(G)
        if G.number_of_nodes() > 34: continue
        h = nx.weisfeiler_lehman_graph_hash(G)
        if h in seen: continue
        seen.add(h); emit('0-núcleo de disco aleatorio' if hard else 'aleatorio', G)
    out.close(); print('discos analizados:', count)
