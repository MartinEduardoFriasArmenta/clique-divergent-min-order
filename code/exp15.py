"""Exhaustive experiment: connected graphs G (not dismantlable) with some H dismantlable and G ->_2 H.
->_2 : H = G - union of private-edge sets of one or more clans R with |R|>2 and eps(R) nonempty.
dismantlable := repeated removal of dominated vertices leaves an edgeless graph.
"""
import sys, time, itertools, subprocess

def parse_g6(line):
    data = [ord(c) - 63 for c in line.strip()]
    n = data[0]; bits = []
    for d in data[1:]:
        for k in range(5, -1, -1): bits.append((d >> k) & 1)
    adj = [0] * n; idx = 0
    for j in range(1, n):
        for i in range(j):
            if bits[idx]: adj[i] |= 1 << j; adj[j] |= 1 << i
            idx += 1
    return adj

def nedges(adj): return sum(bin(a).count('1') for a in adj) // 2

def dismantlable(adj):
    n = len(adj); alive = (1 << n) - 1
    changed = True
    while changed:
        changed = False
        for v in range(n):
            if not (alive >> v) & 1: continue
            Nv = (adj[v] | (1 << v)) & alive
            m = adj[v] & alive
            while m:
                b = m & -m; w = b.bit_length() - 1; m ^= b
                if Nv & ~((adj[w] | b) & alive) == 0:
                    alive &= ~(1 << v); changed = True; break
    for v in range(n):
        if (alive >> v) & 1 and adj[v] & alive: return False
    return True

def cliques(adj):
    n = len(adj); out = []
    def bk(Rm, P, X):
        if P == 0 and X == 0: out.append(Rm); return
        PX = P | X; m = PX; best = -1; piv = 0
        while m:
            b = m & -m; u = b.bit_length() - 1; m ^= b
            c = bin(P & adj[u]).count('1')
            if c > best: best = c; piv = u
        cand = P & ~adj[piv]
        while cand:
            b = cand & -cand; v = b.bit_length() - 1; cand ^= b
            bk(Rm | b, P & adj[v], X & adj[v]); P &= ~b; X |= b
    bk(0, (1 << n) - 1, 0)
    return out

def free_pairs(adj):
    C = cliques(adj); res = []
    for R in C:
        if bin(R).count('1') < 3: continue
        verts = [i for i in range(len(adj)) if (R >> i) & 1]
        eps = []
        for u, v in itertools.combinations(verts, 2):
            m = (1 << u) | (1 << v)
            if not any((Q & m) == m for Q in C if Q != R): eps.append((u, v))
        if eps: res.append((R, eps))
    return res

def remove_edges(adj, eps):
    adj = list(adj)
    for u, v in eps: adj[u] &= ~(1 << v); adj[v] &= ~(1 << u)
    return adj

def analyze(adj, maxsub=14):
    if dismantlable(adj): return None
    fp = free_pairs(adj)
    if not fp: return False
    k_max = len(fp) if len(fp) <= maxsub else 1
    for k in range(1, k_max + 1):
        for sub in itertools.combinations(fp, k):
            eps = [e for _, es in sub for e in es]
            if dismantlable(remove_edges(adj, eps)):
                return [(R, es) for R, es in sub]
    if len(fp) > maxsub:  # try the full union as well
        eps = [e for _, es in fp for e in es]
        if dismantlable(remove_edges(adj, eps)): return list(fp)
    return False

if __name__ == "__main__":
    geng = sys.argv[1]; deadline = time.time() + float(sys.argv[2]) * 60
    out = open(sys.argv[3], 'w'); log = open(sys.argv[4], 'w')
    for n in range(6, 16):
        p = subprocess.Popen([geng, '-c', '-q', str(n)], stdout=subprocess.PIPE, text=True, bufsize=1 << 16)
        tot = nd = pos = 0; t0 = time.time(); complete = True
        stats = {}; written = 0
        for line in p.stdout:
            if time.time() > deadline: complete = False; p.kill(); break
            adj = parse_g6(line); tot += 1
            r = analyze(adj)
            if r is None: continue
            nd += 1
            if r:
                pos += 1
                whole = all(len(es) == bin(R).count('1') * (bin(R).count('1') - 1) // 2 for R, es in r)
                key = ('k=%d' % len(r), 'todas las aristas de R' if whole else 'subconjunto propio')
                stats[key] = stats.get(key, 0) + 1
                if n <= 8 or written < 400 or (not whole and written < 2000):
                    out.write("%d\t%s\t%s\t%s\n" % (n, line.strip(), key[1], repr([(R, es) for R, es in r]))); out.flush(); written += 1
        log.write("n=%d graphs=%d nondismantlable=%d positive=%d time=%.0fs complete=%s stats=%s\n" % (n, tot, nd, pos, time.time() - t0, complete, stats)); log.flush()
        if not complete: break
    out.close(); log.close()
