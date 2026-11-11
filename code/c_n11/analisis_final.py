"""Análisis de los contraejemplos de 11 vértices (estado NO): suspensión o no, 0-núcleo, grado mínimo,
f-vector, chi y números de Betti (Z_2) de Delta(K(G)), y retracción de K(G) a un octaedro O_4 u O_3
(que prueba K-divergencia por el teorema de Neumann-Lara). Escribe contraejemplos_n11.tsv y analisis_n11.tex.
Uso: python3 analisis_final.py <dir_resultados> <dir_salida>
"""
import sys, os, glob, gzip, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import networkx as nx
from homol import betti_z2
from retract import retraction

d, out = sys.argv[1], sys.argv[2]
A = nx.from_graph6_bytes(b'J~Z}eY~qzd_'); L = nx.from_graph6_bytes(b'JrqipkN~~~?')

def issusp(G):
    n = G.number_of_nodes()
    return any(G.degree(u) == n - 2 and G.degree(w) == n - 2 and not G.has_edge(u, w) for u, w in itertools.combinations(G, 2))
def core0(G):
    G = G.copy(); ch = True
    while ch:
        ch = False
        for v in list(G):
            Nv = set(G[v]) | {v}
            if any(Nv <= set(G[w]) | {w} for w in G[v]): G.remove_node(v); ch = True; break
    return G
def clique_graph(G):
    C = [frozenset(c) for c in nx.find_cliques(G)]
    K = nx.Graph(); K.add_nodes_from(range(len(C)))
    K.add_edges_from((a, b) for a, b in itertools.combinations(range(len(C)), 2) if C[a] & C[b]); return K, C
def octahedra(K, k):
    """octaedros inducidos O_k: k parejas no adyacentes, todo lo demás adyacente (búsqueda con retroceso)."""
    V = sorted(K); nonedges = [(u, v) for u, v in itertools.combinations(V, 2) if not K.has_edge(u, v)]
    res = []
    def rec(start, chosen):
        if len(chosen) == 2 * k: res.append(tuple(chosen)); return
        for i in range(start, len(nonedges)):
            u, v = nonedges[i]
            if u in chosen or v in chosen: continue
            if all(K.has_edge(u, x) and K.has_edge(v, x) for x in chosen): rec(i + 1, chosen + [u, v])
    rec(0, []); return res
def retracto(K):
    for k in (4, 3):
        for S in octahedra(K, k):
            if retraction(K, S): return k
    return None

rows = []
for f in sorted(glob.glob(os.path.join(d, 'parte_*.tsv.gz'))):
    for line in gzip.open(f, 'rt'):
        g6, m, st = line.rstrip('\n').split('\t')
        if st == 'NO': rows.append(g6)
import ast, re as _re
cache = {}   # resultados previos (archivo de bitácora de una corrida anterior), para reanudar
if len(sys.argv) > 3 and os.path.exists(sys.argv[3]):
    for line in open(sys.argv[3]):
        m = _re.match(r'^(\S+) (\d+) (susp|NO-susp) (\d+) (-?\d+) (\[.*\]) ?(O_\d)? ?(\S*)$', line.rstrip('\n'))
        if m: cache[m.group(1)] = (int(m.group(4)), int(m.group(5)), ast.literal_eval(m.group(6)), m.group(7) or '')
res = []
for g6 in rows:
    G = nx.from_graph6_bytes(g6.encode())
    if g6 in cache:
        nK, chi, b, rt = cache[g6]
        hom = ', '.join('$\\tilde H_{%d}=\\mathbb Z_2^{%d}$' % (i, x) if x > 1 else '$\\tilde H_{%d}=\\mathbb Z_2$' % i for i, x in enumerate(b) if x) or '0'
        com = 'A' if nx.is_isomorphic(G, A) else ('G_LPV' if nx.is_isomorphic(G, L) else '')
        res.append(dict(g6=g6, aristas=G.number_of_edges(), mindeg=min(dd for _, dd in G.degree()), susp=issusp(G), nucleo=core0(G).number_of_nodes(),
                        nK=nK, chi=chi, betti=b, hom=hom, retracto=rt, com=com))
        continue
    K, C = clique_graph(G)
    fv = {}
    for c in nx.enumerate_all_cliques(K): fv[len(c)] = fv.get(len(c), 0) + 1
    chi = sum((-1) ** (k - 1) * v for k, v in fv.items())
    b = betti_z2(K, maxdim=max(fv) + 1)[:max(fv)]
    hom = ', '.join('$\\tilde H_{%d}=\\mathbb Z_2^{%d}$' % (i, x) if x > 1 else '$\\tilde H_{%d}=\\mathbb Z_2$' % i for i, x in enumerate(b) if x) or '0'
    ok = retracto(K)
    com = 'A' if nx.is_isomorphic(G, A) else ('G_LPV' if nx.is_isomorphic(G, L) else '')
    r = dict(g6=g6, aristas=G.number_of_edges(), mindeg=min(dd for _, dd in G.degree()), susp=issusp(G),
             nucleo=core0(G).number_of_nodes(), nK=K.number_of_nodes(), chi=chi, betti=b, hom=hom,
             retracto=('O_%d' % ok) if ok else '', com=com)
    res.append(r); print(r['g6'], r['aristas'], 'susp' if r['susp'] else 'NO-susp', r['nK'], r['chi'], r['betti'], r['retracto'], com, flush=True)

res.sort(key=lambda r: (r['susp'], r['aristas'], r['g6']))
with open(os.path.join(out, 'contraejemplos_n11.tsv'), 'w') as fo:
    fo.write('graph6\taristas\tgrado_min\tsuspension\tnucleo0\t|K|\tchi(K)\tbetti_Z2_reducidos(K)\tretracto_octaedrico\tcomentario\n')
    for r in res:
        fo.write('%s\t%d\t%d\t%s\t%d\t%d\t%d\t%s\t%s\t%s\n' % (r['g6'], r['aristas'], r['mindeg'], 'si' if r['susp'] else 'no', r['nucleo'], r['nK'], r['chi'], r['betti'], r['retracto'] or '-', r['com']))
def esc(s): return s.replace('\\', '\\textbackslash{}').replace('_', '\\_').replace('~', '\\textasciitilde{}').replace('^', '\\textasciicircum{}').replace('{', '\\{').replace('}', '\\}').replace('#', '\\#').replace('&', '\\&').replace('%', '\\%')
ns = [r for r in res if not r['susp']]; su = [r for r in res if r['susp']]
T = [r'\section*{Análisis de los contraejemplos}',
     'Hay %d gráficas de $11$ vértices que son $1$-desmantelables y no $0$-desmantelables con $K(G)$ no $1$-desmantelable. De ellas, %d son suspensiones (las de la Proposición 5.2) y \\textbf{%d no son suspensiones}, lo que responde afirmativamente la Pregunta 7.1. ' % (len(res), len(su), len(ns)) + ('Ninguna tiene vértices $0$-dominados. ' if all(r['nucleo'] == 11 for r in res) else '%d de ellas no tienen vértices $0$-dominados. ' % sum(1 for r in res if r['nucleo'] == 11)),
     'En %d de las %d, $K(G)$ se retrae a un octaedro ($O_4$ u $O_3$), de modo que $G$ es $K$-divergente por el teorema de Neumann-Lara; entre las no suspensiones, %d de %d. ' % (sum(1 for r in res if r['retracto']), len(res), sum(1 for r in ns if r['retracto']), len(ns)),
     'La característica de Euler de $\\Delta(K(G))$ es distinta de $1$ en %d casos, lo que da una prueba independiente (Lema 2.2) de que $K(G)$ no es $1$-desmantelable.' % sum(1 for r in res if r['chi'] != 1),
     r'''
\begin{longtable}{rlrrrrll}\toprule
\# & graph6 & aristas & grado mín. & $|K(G)|$ & $\chi(\Delta K)$ & homología reducida de $\Delta(K(G))$ ($\mathbb Z_2$) & retracto\\\midrule\endhead
\multicolumn{8}{l}{\emph{No suspensiones}}\\''']
for i, r in enumerate(ns, 1):
    T.append('%d & \\texttt{%s} & %d & %d & %d & %d & %s & %s\\\\' % (i, esc(r['g6']), r['aristas'], r['mindeg'], r['nK'], r['chi'], r['hom'], '$%s$' % r['retracto'] if r['retracto'] else '—'))
T.append(r'\midrule\multicolumn{8}{l}{\emph{Suspensiones (Proposición 5.2)}}\\')
for i, r in enumerate(su, len(ns) + 1):
    T.append('%d & \\texttt{%s}%s & %d & %d & %d & %d & %s & %s\\\\' % (i, esc(r['g6']), (' ($%s$)' % r['com'].replace('G_LPV', 'G_{LPV}')) if r['com'] else '', r['aristas'], r['mindeg'], r['nK'], r['chi'], r['hom'], '$%s$' % r['retracto'] if r['retracto'] else '—'))
T.append(r'\bottomrule\end{longtable}')
open(os.path.join(out, 'analisis_n11.tex'), 'w').write('\n'.join(T) + '\n')
print('suspensiones', len(su), 'no suspensiones', len(ns))
