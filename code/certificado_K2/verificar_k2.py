"""Verificación independiente del certificado: K^2(G) se retrae a un octaedro O_16, G = JEhfrz{~~~?."""
import json, itertools
import networkx as nx
d = json.load(open('k2_retraccion_JEhfrz.json'))
G = nx.from_graph6_bytes(b'JEhfrz{~~~?')
# reconstruir K(G) y K^2(G) desde cero y compararlos con los clanes guardados
CX = sorted(sorted(c) for c in nx.find_cliques(G))
assert CX == sorted(d['clanes_KG']), 'clanes de G distintos'
CX = [frozenset(c) for c in d['clanes_KG']]
X = nx.Graph(); X.add_nodes_from(range(len(CX))); X.add_edges_from((a, b) for a, b in itertools.combinations(range(len(CX)), 2) if CX[a] & CX[b])
CY2 = sorted(sorted(c) for c in nx.find_cliques(X)); assert CY2 == sorted(d['clanes_K2G']), 'clanes de K(G) distintos'
CY = [frozenset(c) for c in d['clanes_K2G']]
Y = nx.Graph(); Y.add_nodes_from(range(len(CY))); Y.add_edges_from((a, b) for a, b in itertools.combinations(range(len(CY)), 2) if CY[a] & CY[b])
S = d['O16']; r = {int(k): v for k, v in d['r'].items()}
print('K(G): %d vértices; K^2(G): %d vértices, %d aristas' % (X.number_of_nodes(), Y.number_of_nodes(), Y.number_of_edges()))
# 1. S induce O_16 con parejas antipodales (S[0],S[1]), (S[2],S[3]), ...
anti = {S[i]: S[i ^ 1] for i in range(32)}
assert len(set(S)) == 32
assert all((not Y.has_edge(u, v)) == (anti[u] == v) for u, v in itertools.combinations(S, 2)), 'S no induce O_16'
# 2. r es identidad en S y está definida en todo K^2(G) con valores en S
assert set(r) == set(Y.nodes()) and all(r[s] == s for s in S) and all(v in S for v in r.values())
# 3. condición de morfismo
bad = [(u, v) for u, v in Y.edges() if r[u] != r[v] and not Y.has_edge(r[u], r[v])]
print('S induce O_16: sí; r es identidad en S: sí; aristas que violan la condición de morfismo:', len(bad))
print('Clanes de K(G) que forman el O_16 (cada clan de K(G) es un conjunto de clanes de G):')
for i in range(0, 32, 2): print('  ', sorted(CY[S[i]]), ' <-> ', sorted(CY[S[i+1]]))
