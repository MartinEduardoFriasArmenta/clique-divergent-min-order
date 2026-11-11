"""Todas las suspensiones Susp(H) con |H| = 9 (geng): ¿cuáles son 1-desmantelables con K(Susp H) no 1-desmantelable?
Uso: python3 susp11.py <geng> <salida.tsv>
"""
import sys, time, subprocess
import networkx as nx
from exp15 import parse_g6, cliques
from discoid import zero_dismantlable, Discoid
geng = sys.argv[1]; out = open(sys.argv[2], 'w')
def kgraph(adj):
    C = cliques(adj); m = len(C); K = [0] * m
    for i in range(m):
        for j in range(i + 1, m):
            if C[i] & C[j]: K[i] |= 1 << j; K[j] |= 1 << i
    return K
p = subprocess.Popen([geng, '-c', '-q', '9'], stdout=subprocess.PIPE, text=True, bufsize=1 << 16)
t0 = time.time(); tot = disc = found = 0; seen = set()
for line in p.stdout:
    tot += 1; adjH = parse_g6(line); SH = (1 << 9) - 1
    if zero_dismantlable(adjH, SH): continue            # Susp(H) sería 0-desmantelable
    if not Discoid(adjH, cap=100000).test(SH): continue  # Susp(H) 1-desm <=> H 1-desm
    disc += 1
    adj = [a | (1 << 9) | (1 << 10) for a in adjH] + [SH, SH]   # ápices 9 y 10
    S = (1 << 11) - 1
    assert Discoid(adj, cap=100000).test(S) and not zero_dismantlable(adj, S)
    K = kgraph(adj); m = len(K); SK = (1 << m) - 1
    r = Discoid(K, cap=600000).test(SK)
    if r is False:
        found += 1
        G = nx.Graph(); G.add_nodes_from(range(11))
        for i in range(11):
            for j in range(i + 1, 11):
                if (adj[i] >> j) & 1: G.add_edge(i, j)
        g6 = nx.to_graph6_bytes(G, header=False).strip().decode()
        out.write('%s\t%s\t%d\t%d\n' % (line.strip(), g6, G.number_of_edges(), m)); out.flush()
        print('H=%s  Susp(H)=%s  aristas=%d  |K|=%d' % (line.strip(), g6, G.number_of_edges(), m), flush=True)
    elif r is None: print('indecidido', line.strip(), flush=True)
print('H conexas=%d, H 1-desm no 0-desm=%d, Susp(H) con K no 1-desm=%d, tiempo=%.0fs' % (tot, disc, found, time.time() - t0))
