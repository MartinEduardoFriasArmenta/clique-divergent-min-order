"""K-nulidad de las gráficas registradas por exh_disc2.py (1-desmantelables, no 0-desmantelables).
Uso: python3 kbehav.py <entrada.tsv> [max_iter=12] [max_orden=200]
Itera K hasta llegar a K_1 (K-nula) o exceder los límites; como el borrado de vértices 0-dominados
preserva el K-comportamiento (Frías-Armenta, Neumann-Lara, Pizaña 2004), en cada paso se sustituye
la gráfica por su 0-núcleo, lo que mantiene pequeños los órdenes.
"""
import sys, time
from exp15 import parse_g6, cliques
from discoid import zero_dominated, popcount
def core0(adj):
    n = len(adj); S = (1 << n) - 1; ch = True
    while ch:
        ch = False; m = S
        while m:
            b = m & -m; v = b.bit_length() - 1; m ^= b
            if zero_dominated(adj, S, v): S &= ~b; ch = True; break
    idx = [v for v in range(n) if (S >> v) & 1]; pos = {v: i for i, v in enumerate(idx)}
    new = [0] * len(idx)
    for i, v in enumerate(idx):
        m = adj[v] & S
        while m:
            b = m & -m; w = b.bit_length() - 1; m ^= b; new[i] |= 1 << pos[w]
    return new
def kgraph(adj):
    C = cliques(adj); m = len(C); K = [0] * m
    for i in range(m):
        for j in range(i + 1, m):
            if C[i] & C[j]: K[i] |= 1 << j; K[j] |= 1 << i
    return K
maxit = int(sys.argv[2]) if len(sys.argv) > 2 else 12; maxn = int(sys.argv[3]) if len(sys.argv) > 3 else 200
t0 = time.time(); stats = {}; hist = {}; unresolved = []
for line in open(sys.argv[1]):
    g6 = line.split('\t')[0]; adj = core0(parse_g6(g6)); status = None
    for it in range(maxit + 1):
        if len(adj) == 1: status = 'K-nula (%d pasos)' % it; hist[it] = hist.get(it, 0) + 1; break
        if len(adj) > maxn: status = 'orden > %d' % maxn; break
        adj = core0(kgraph(adj))
    if status is None: status = 'sin decidir en %d pasos' % maxit
    key = status.split(' (')[0]; stats[key] = stats.get(key, 0) + 1
    if not status.startswith('K-nula'): unresolved.append((g6, status))
print(sys.argv[1], stats, 'pasos hasta K_1 (tras 0-núcleos):', dict(sorted(hist.items())), 'tiempo=%.0fs' % (time.time() - t0))
for g6, st in unresolved[:50]: print('  ', g6, st)
