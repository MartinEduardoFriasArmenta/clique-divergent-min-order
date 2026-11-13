"""Verificación de las propiedades de la gráfica A (y del ejemplo de Larrión-Pizaña-Villarroel-Flores).
Uso: python3 propiedades_A.py [graph6]      (por omisión A = J~Z}eY~qzd_)
Comprueba:
  1. A no es 0-desmantelable, sí es 1-desmantelable (orden explícito).
  2. A es strong vertex I-contractible (Algoritmo 1 de Dochtermann-Espinoza-Frías-Hernández 2023, y versión exacta).
  3. Delta(A) es colapsible: sucesión explícita de colapsos elementales hasta un vértice.
  4. chi(Delta(A)) = 1 y homología reducida (Z_2) nula.
  5. K(A): orden, tamaño, no 0-desmantelable, NO 1-desmantelable (búsqueda exacta), chi(Delta(K(A))) = -1, Betti Z_2.
  6. K(A) se retrae a un octaedro O_4  =>  A es K-divergente (Neumann-Lara).
  7. A = Susp(H) con H = A - {u,w} para dos vértices no adyacentes u,w de grado n-2 (si existen).
"""
import sys, itertools, random
import networkx as nx
from functools import lru_cache
from discoid import to_adj, popcount, zero_dominated, zero_dismantlable, one_dominated, Discoid, is_discoid, clique_graph
from homol import betti_z2
from retract import octahedra, retraction

g6 = sys.argv[1] if len(sys.argv) > 1 else 'J~Z}eY~qzd_'
G = nx.from_graph6_bytes(g6.encode()); n = G.number_of_nodes()
print('Gráfica %s: %d vértices, %d aristas' % (g6, n, G.number_of_edges()))
adj, nodes = to_adj(G); S = (1 << n) - 1

# 1. desmantelabilidad
print('1. 0-desmantelable:', zero_dismantlable(adj, S))
ok, order = is_discoid(G, with_order=True)
print('   1-desmantelable:', ok, ' orden de borrado:', order)
if ok:  # verificación independiente del orden
    T = S
    for v in order:
        i = nodes.index(v); assert one_dominated(adj, T, i), v; T &= ~(1 << i)
    assert popcount(T) == 1
    print('   (orden verificado: cada vértice es 1-dominado en el momento de borrarlo)')

# 2. strong vertex I-contractible
def SIV_alg1(H):
    """Algoritmo 1 de DEFH 2023 (voraz: toma el primer vértice I_SV-contráctil)."""
    if H.number_of_nodes() == 0: return False
    if H.number_of_nodes() == 1: return True
    for v in list(H):
        if SIV_alg1(H.subgraph(H[v])): return SIV_alg1(H.subgraph([u for u in H if u != v]))
    return False
@lru_cache(maxsize=None)
def SIV_exact(mask):
    """G[mask] en I_SV (exacto: existe algún orden)."""
    if popcount(mask) == 1: return True
    m = mask
    while m:
        b = m & -m; v = b.bit_length() - 1; m ^= b
        N = adj[v] & mask
        if N and SIV_exact(N) and SIV_exact(mask & ~b): return True
    return False
print('2. strong vertex I-contractible (Algoritmo 1, voraz):', SIV_alg1(G), '; (exacto):', SIV_exact(S))

# 3. colapsibilidad de Delta(G)
def collapse_sequence(G, rnd, tries=50):
    simp = set(frozenset(c) for c in nx.enumerate_all_cliques(G))
    for t in range(tries):
        cx = set(simp); seq = []
        while len(cx) > 1:
            maximal = [s for s in cx if not any(s < u for u in cx)]
            free = []
            for tau in maximal:
                if len(tau) == 1: continue
                for v in tau:
                    sigma = tau - {v}
                    if not any(sigma < u and u != tau for u in cx): free.append((sigma, tau))
            if not free: break
            sigma, tau = rnd.choice(free); cx.discard(sigma); cx.discard(tau); seq.append((sorted(sigma), sorted(tau)))
        if len(cx) == 1: return seq
    return None
rnd = random.Random(0); seq = collapse_sequence(G, rnd)
fvec = {}
for c in nx.enumerate_all_cliques(G): fvec[len(c)] = fvec.get(len(c), 0) + 1
print('3. f-vector de Delta(G):', [fvec[k] for k in sorted(fvec)], ' colapsible:', seq is not None, '(%d colapsos elementales)' % (len(seq) if seq else 0))
if seq:
    print('   primeros colapsos (cara libre -> cara maximal):', ', '.join('%s<%s' % (''.join(map(str, a)), ''.join(map(str, b))) for a, b in seq[:8]), '...')
    open('colapsos_%s.txt' % g6.replace('~', 'T').replace('\\', 'B').replace('}', 'R').replace('_', 'U').replace('?', 'Q').replace('{', 'L').replace('|', 'P').replace('^', 'C').replace('`', 'G'), 'w').write('\n'.join('%s %s' % (' '.join(map(str, a)), ' '.join(map(str, b))) for a, b in seq) + '\n')

# 4. homología
chi = sum((-1) ** (k - 1) * v for k, v in fvec.items())
print('4. chi(Delta(G)) =', chi, '; Betti reducidos Z_2 =', betti_z2(G, maxdim=max(fvec) + 1)[:max(fvec)])

# 5. K(G)
K, C = clique_graph(G); m = K.number_of_nodes(); adjK, _ = to_adj(K); SK = (1 << m) - 1
print('5. K(G): %d vértices (clanes), %d aristas; clanes:' % (m, K.number_of_edges()), [''.join(map(str, sorted(c))) for c in C])
print('   K(G) 0-desmantelable:', zero_dismantlable(adjK, SK))
DK = Discoid(adjK, cap=3000000); rK = DK.test(SK)
print('   K(G) 1-desmantelable:', rK, '(búsqueda exacta, %d subconjuntos explorados)' % len(DK.memo))
fK = {}
for c in nx.enumerate_all_cliques(K): fK[len(c)] = fK.get(len(c), 0) + 1
chiK = sum((-1) ** (k - 1) * v for k, v in fK.items())
print('   f-vector de Delta(K(G)):', [fK[k] for k in sorted(fK)], ' chi =', chiK, '; Betti reducidos Z_2 =', betti_z2(K, maxdim=max(fK) + 1)[:max(fK)])
d1K = [i for i in range(m) if one_dominated(adjK, SK, i)]
print('   vértices 1-dominados de K(G): %d de %d' % (len(d1K), m))

# 6. retracto octaédrico
found = None
for k in (4, 3):
    for Sx in octahedra(K, k):
        r = retraction(K, Sx)
        if r: found = (k, Sx, r); break
    if found: break
if found:
    k, Sx, r = found
    print('6. K(G) se retrae al octaedro O_%d con vértices' % k, [''.join(map(str, sorted(C[i]))) for i in Sx])
    print('   retracción (clan -> imagen):', {''.join(map(str, sorted(C[v]))): ''.join(map(str, sorted(C[x]))) for v, x in r.items() if v != x})
    bad = [(u, v) for u, v in K.edges() if r[u] != r[v] and not K.has_edge(r[u], r[v])]
    print('   verificación de la condición de morfismo en las %d aristas de K(G): violaciones = %d' % (K.number_of_edges(), len(bad)))
    print('   => G es K-divergente (Neumann-Lara 1978).')
else:
    print('6. no se encontró retracto octaédrico en K(G)')

# 7. estructura de suspensión
apex = [v for v in G if G.degree(v) == n - 2]
pairs = [(u, w) for u, w in itertools.combinations(apex, 2) if not G.has_edge(u, w)]
if pairs:
    u, w = pairs[0]; H = G.subgraph([v for v in G if v not in (u, w)])
    print('7. G = Susp(H) con ápices %s,%s; H: %d vértices, %d aristas, plana=%s, graph6=%s' % (u, w, H.number_of_nodes(), H.number_of_edges(), nx.check_planarity(H)[0], nx.to_graph6_bytes(H, header=False).strip().decode()))
    adjH, _ = to_adj(H); SH = (1 << H.number_of_nodes()) - 1
    print('   H 0-desmantelable:', zero_dismantlable(adjH, SH), '; H 1-desmantelable:', is_discoid(H))
else:
    print('7. G no es una suspensión')
