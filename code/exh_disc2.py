"""Enumeración exhaustiva (geng): gráficas conexas con n vértices que son 1-desmantelables y no 0-desmantelables.
Para cada una se comprueba si K(G) es 1-desmantelable y se guarda su graph6 en <salida>.
Uso: python3 exh_disc2.py <geng> <n> <salida> [res/mod]
"""
import sys, time, subprocess
from exp15 import parse_g6, cliques
from discoid import zero_dismantlable, Discoid
geng = sys.argv[1]; n = int(sys.argv[2]); out = open(sys.argv[3], 'w'); extra = sys.argv[4:]
def kgraph(adj):
    C = cliques(adj); m = len(C); K = [0] * m
    for i in range(m):
        for j in range(i + 1, m):
            if C[i] & C[j]: K[i] |= 1 << j; K[j] |= 1 << i
    return K
p = subprocess.Popen([geng, '-c', '-q', str(n)] + extra, stdout=subprocess.PIPE, text=True, bufsize=1 << 16)
t0 = time.time(); tot = disc = kndz = bad = und = 0
for line in p.stdout:
    tot += 1; adj = parse_g6(line); S = (1 << n) - 1
    if zero_dismantlable(adj, S): continue
    if not Discoid(adj, cap=100000).test(S): continue
    disc += 1
    K = kgraph(adj); m = len(K); SK = (1 << m) - 1
    if zero_dismantlable(K, SK): kd = '0'
    else:
        kndz += 1; r = Discoid(K, cap=400000).test(SK)
        kd = '1' if r else ('NO' if r is False else '?')
        if r is False: bad += 1
        if r is None: und += 1
    out.write('%s\t%d\t%s\n' % (line.strip(), m, kd)); out.flush()
print('n=%d conexas=%d 1-desm_no_0-desm=%d K_no_0-desm=%d K_no_1-desm=%d indecidido=%d tiempo=%.0fs' % (n, tot, disc, kndz, bad, und, time.time() - t0), flush=True)
