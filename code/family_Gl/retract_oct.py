"""Retracción de una gráfica X a un octaedro inducido S (lista de vértices en parejas antipodales s0,s1, s2,s3, ...)
como problema de satisfacción de restricciones: r(u) y r(v) no pueden ser antipodales si uv es arista.
Búsqueda con dominio mínimo primero y verificación hacia adelante. Devuelve el diccionario r o None."""
import sys
sys.setrecursionlimit(10000)
def retract_octahedron(X, S, limit=2_000_000):
    S = list(S); anti = {}
    for i in range(0, len(S), 2): anti[S[i]] = S[i+1]; anti[S[i+1]] = S[i]
    r = {s: s for s in S}
    dom = {}
    for x in X:
        if x in r: continue
        bad = {anti[s] for s in X[x] if s in r}
        dom[x] = set(S) - bad
        if not dom[x]: return None
    nodes = [0]
    def assign(x, val, doms):
        new = {}
        for y in X[x]:
            if y in doms:
                if anti[val] in doms[y]:
                    d = doms[y] - {anti[val]}
                    if not d: return None
                    new[y] = d
        return new
    def bt(doms):
        nodes[0] += 1
        if nodes[0] > limit: raise TimeoutError
        if not doms: return True
        x = min(doms, key=lambda v: (len(doms[v]), -len(X[v])))
        for val in sorted(doms[x]):
            upd = assign(x, val, doms)
            if upd is None: continue
            nd = {y: (upd[y] if y in upd else d) for y, d in doms.items() if y != x}
            r[x] = val
            if bt(nd): return True
            del r[x]
        return False
    try:
        return dict(r) if bt(dom) else None
    except TimeoutError:
        return 'limite'

def retract_octahedron_fixed(X, S, fixed, limit=2_000_000):
    """Como retract_octahedron, pero con valores prefijados fixed[v] (en S) para algunos vértices v fuera de S."""
    S = list(S); anti = {}
    for i in range(0, len(S), 2): anti[S[i]] = S[i+1]; anti[S[i+1]] = S[i]
    r = {s: s for s in S}; r.update(fixed)
    for u, v in X.edges():                       # consistencia de lo prefijado
        if u in r and v in r and r[u] == anti[r[v]]: return None
    dom = {}
    for x in X:
        if x in r: continue
        bad = {anti[r[s]] for s in X[x] if s in r}
        dom[x] = set(S) - bad
        if not dom[x]: return None
    nodes = [0]
    def bt(doms):
        nodes[0] += 1
        if nodes[0] > limit: raise TimeoutError
        if not doms: return True
        x = min(doms, key=lambda v: (len(doms[v]), -len(X[v])))
        for val in sorted(doms[x]):
            nd = {}; okk = True
            for y, d in doms.items():
                if y == x: continue
                if y in X[x] and anti[val] in d:
                    d = d - {anti[val]}
                    if not d: okk = False; break
                nd[y] = d
            if not okk: continue
            r[x] = val
            if bt(nd): return True
            del r[x]
        return False
    try:
        return dict(r) if bt(dom) else None
    except TimeoutError:
        return 'limite'
