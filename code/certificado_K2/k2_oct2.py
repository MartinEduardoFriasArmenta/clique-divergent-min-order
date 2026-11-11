import itertools, sys, time, json
import networkx as nx
sys.path.insert(0, '.')
from retract_oct import retract_octahedron
exec(open('k2_oct.py').read().split("G = nx.from_graph6_bytes")[0].split("from retract import retraction")[1])
G = nx.from_graph6_bytes(b'JEhfrz{~~~?')
X, CX = K(G); Y, CY = K(X)
# control: el método encuentra los retractos O_5 conocidos de otro ejemplo y ninguno O_5 en K(G) aquí
O5 = octahedra(X, 5)
print('O_5 de K(G) que son retracto (control, se espera 0):', sum(1 for S in O5 if retract_octahedron(X, S) not in (None,)), flush=True)
t0 = time.time(); found = 0; tested = 0
for S in O5:
    pairs = [(S[2*i], S[2*i+1]) for i in range(5)]
    bitsl = list(itertools.product((0, 1), repeat=5))
    choices = [[i for i, c in enumerate(CY) if frozenset(p[b] for p, b in zip(pairs, bits)) <= c] for bits in bitsl]
    for combo in itertools.product(*choices):
        ok = True
        for a in range(32):
            for b in range(a + 1, 32):
                anti = all(x != y for x, y in zip(bitsl[a], bitsl[b]))
                if anti == bool(CY[combo[a]] & CY[combo[b]]): ok = False; break
            if not ok: break
        if not ok: continue
        found += 1
        order = []
        for a in range(32):
            b = bitsl.index(tuple(1 - x for x in bitsl[a]))
            if a < b: order += [combo[a], combo[b]]
        r = retract_octahedron(Y, order); tested += 1
        if r == 'limite': print('límite alcanzado', flush=True); continue
        if r:
            bad = [(u, v) for u, v in Y.edges() if r[u] != r[v] and not Y.has_edge(r[u], r[v])]
            print('RETRACCIÓN K^2(G) -> O_16 encontrada; violaciones =', len(bad), round(time.time() - t0), 's', flush=True)
            json.dump({'O16': order, 'r': {str(k): v for k, v in r.items()}, 'clanes_KG': [sorted(c) for c in CX], 'clanes_K2G': [sorted(c) for c in CY]}, open('k2_retraccion_JEhfrz.json', 'w'))
            sys.exit(0)
print('O_16 inducidos:', found, 'probados:', tested, '; ninguno es retracto', round(time.time() - t0), 's')
