/* gen_disks_ca.c -- Whitney triangulations of the disk (flag disks) on exactly N vertices, up to isomorphism,
 * by McKay's canonical construction path.  Unlike gen_disks.c it stores nothing, so it reaches larger N, and the
 * work can be split into independent parts (res/mod), as in geng.
 *
 * Construction (see gen_disks.c): a child G is obtained from a parent P by adding a vertex v adjacent to an induced
 * path x_0...x_k (k >= 1) of consecutive boundary vertices of P.  A boundary vertex u of G is "removable" if no
 * interior vertex of the path N(u) lies on the boundary of G; then G - u is a flag disk (Proposition 2.7 of the paper),
 * and every flag disk with >= 4 vertices has a removable vertex.  The canonical removable vertex m(G) is chosen with
 * an invariant (degree, sum of the degrees of the neighbors) and, among the removable vertices with maximal
 * invariant, by nauty's canonical labelling.  G is accepted iff v is in the orbit of m(G) under Aut(G); isomorphic
 * children of the same parent are removed with a small per-parent table.
 *
 * Usage:  gen_disks_ca N [res mod [splitlevel]]     (graph6 on stdout; default splitlevel = min(N-1, 11))
 *         writes the flag disks on N vertices; with res/mod only those descending from the parents of order
 *         splitlevel whose index is congruent to res modulo mod.  Counts are printed to stderr.
 * Build:  gcc -O3 -march=native -o gen_disks_ca gen_disks_ca.c -DWORDSIZE=64 -DMAXN=WORDSIZE -lnautyL1
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nauty.h"
#include "gtools.h"

static int N, RES = 0, MOD = 1, SPLIT;
static long long splitcount = 0, outcount = 0;

static int boundary(const graph *g, int n, int *cyc, int *onb) {
    setword bnb[MAXN];
    for (int i = 0; i < n; i++) { bnb[i] = 0; onb[i] = 0; }
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (ISELEMENT(&g[i], j) && POPCOUNT(g[i] & g[j]) == 1) { bnb[i] |= bit[j]; bnb[j] |= bit[i]; onb[i] = onb[j] = 1; }
    int start = -1; for (int i = 0; i < n; i++) if (onb[i]) { start = i; break; }
    int b = 0, prev = -1, cur = start;
    do { cyc[b++] = cur; int nxt = -1;
         for (int j = 0; j < n; j++) if (ISELEMENT(&bnb[cur], j) && j != prev) { nxt = j; break; }
         prev = cur; cur = nxt; } while (cur != start && cur >= 0 && b <= n);
    return b;
}

static void rec(graph *g, int n) {
    if (n == N) { fputs(ntog6(g, 1, n), stdout); outcount++; return; }
    if (n == SPLIT) { if ((splitcount++ % MOD) != RES) return; }
    int cyc[MAXN], onb[MAXN];
    int b = boundary(g, n, cyc, onb);
    graph h[MAXN], cg[MAXN];
    graph *kids = NULL; int nk = 0, capk = 0;                 /* canonical forms of accepted children */
    int lab[MAXN], ptn[MAXN], orbits[MAXN], deg[MAXN], inv[MAXN], rem[MAXN], hb[MAXN];
    DEFAULTOPTIONS_GRAPH(options); statsblk stats;
    options.getcanon = TRUE; options.defaultptn = FALSE;
    int m = n + 1, v = n;
    for (int s = 0; s < b; s++)
        for (int k = 1; k <= b - 1; k++) {
            int ok = 1;
            for (int i = 0; i <= k && ok; i++)
                for (int j = i + 2; j <= k; j++)
                    if (ISELEMENT(&g[cyc[(s + i) % b]], cyc[(s + j) % b])) { ok = 0; break; }
            if (!ok) break;
            for (int i = 0; i < n; i++) h[i] = g[i];
            h[v] = 0;
            for (int i = 0; i <= k; i++) { int x = cyc[(s + i) % b]; h[x] |= bit[v]; h[v] |= bit[x]; }
            /* boundary of h: that of g with x_1..x_{k-1} replaced by v */
            for (int i = 0; i < n; i++) hb[i] = onb[i];
            for (int i = 1; i < k; i++) hb[cyc[(s + i) % b]] = 0;
            hb[v] = 1;
            for (int u = 0; u < m; u++) deg[u] = POPCOUNT(h[u]);
            int best = -1;
            for (int u = 0; u < m; u++) {
                rem[u] = 0;
                if (!hb[u]) continue;
                /* removable: the neighbors of u other than its two boundary neighbors are interior */
                int nb = 0; setword nu = h[u];
                while (nu) { int w = FIRSTBITNZ(nu); nu ^= bit[w]; if (hb[w]) nb++; }
                /* u's boundary neighbors are exactly its neighbors joined to u by a boundary edge */
                int bn = 0; nu = h[u];
                while (nu) { int w = FIRSTBITNZ(nu); nu ^= bit[w]; if (POPCOUNT(h[u] & h[w]) == 1) bn++; }
                if (nb != bn) continue;                       /* some interior path vertex lies on the boundary */
                rem[u] = 1;
                int sd = 0; nu = h[u]; while (nu) { int w = FIRSTBITNZ(nu); nu ^= bit[w]; sd += deg[w]; }
                inv[u] = (deg[u] << 12) | sd;
                if (inv[u] > best) best = inv[u];
            }
            if (!rem[v] || inv[v] != best) continue;           /* v cannot be the canonical removable vertex */
            int cnt = 0, single = 1;
            for (int u = 0; u < m; u++) if (rem[u] && inv[u] == best) { cnt++; if (u != v) single = 0; }
            /* colouring: removable vertices of maximal invariant first, then the rest (intrinsic, hence canonical) */
            int p = 0;
            for (int u = 0; u < m; u++) if (rem[u] && inv[u] == best) lab[p++] = u;
            int c1 = p;
            for (int u = 0; u < m; u++) if (!(rem[u] && inv[u] == best)) lab[p++] = u;
            for (int i = 0; i < m; i++) ptn[i] = 1;
            ptn[c1 - 1] = 0; ptn[m - 1] = 0;
            densenauty(h, lab, ptn, orbits, &options, &stats, 1, m, cg);
            if (!single) {
                int mg = lab[c1 - 1];                          /* canonical removable vertex: last of its cell */
                if (orbits[v] != orbits[mg]) continue;
            }
            int dup = 0;
            for (int i = 0; i < nk && !dup; i++) if (!memcmp(kids + (size_t)i * m, cg, m * sizeof(graph))) dup = 1;
            if (dup) continue;
            if (nk == capk) { capk = capk ? 2 * capk : 32; kids = realloc(kids, (size_t)capk * m * sizeof(graph)); }
            memcpy(kids + (size_t)nk * m, cg, m * sizeof(graph)); nk++;
            rec(h, m);
        }
    free(kids);
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: gen_disks_ca N [res mod [splitlevel]]\n"); return 1; }
    N = atoi(argv[1]);
    if (argc >= 4) { RES = atoi(argv[2]); MOD = atoi(argv[3]); }
    SPLIT = argc >= 5 ? atoi(argv[4]) : (N - 1 < 11 ? N - 1 : 11);
    if (N < 3 || N > MAXN - 1) { fprintf(stderr, "N out of range\n"); return 1; }
    nauty_check(WORDSIZE, 1, MAXN, NAUTYVERSIONID);
    graph t[3] = {bit[1] | bit[2], bit[0] | bit[2], bit[0] | bit[1]};
    if (MOD > 1 && SPLIT < 3) SPLIT = 3;
    rec(t, 3);
    fprintf(stderr, ">Z %lld flag disks on %d vertices (part %d/%d)\n", outcount, N, RES, MOD);
    return 0;
}
