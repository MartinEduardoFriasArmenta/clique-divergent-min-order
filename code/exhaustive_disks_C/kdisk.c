/* kdisk.c -- iterated clique graphs of Whitney disks: 0- and 1-dismantlability of the iterates (Section 7).
 *
 * For each input graph B (graph6, one per line) it computes H_0 = 0-core(B) and H_i = 0-core(K(H_{i-1})),
 * i = 1..IT.  By Frias-Armenta--Neumann-Lara--Pizana (2004), K^i(B) 0-dismantles to H_i, hence
 *   - K^i(B) is 0-dismantlable  <=>  H_i = K_1   (and then B is K-null);
 *   - if H_i is 1-dismantlable, so is K^i(B)     (delete first the dominated vertices, which are 1-dominated).
 * Working with the cores instead of the full iterates keeps the graphs small.
 *
 * Output (tab separated):  graph6  n  m  e_0 e_1 ... [K1 | >L]
 *   e_0 = "n0:c0:s0" with n0 = |B|, c0 = |H_0|; e_i = "ni:ci:si" with ni = |K(H_{i-1})|, ci = |H_i|;
 *   s = 0  (H_i = K_1, i.e. K^i(B) is 0-dismantlable),  1 (H_i is 1-dismantlable, certificate by greedy or exact search),
 *       NO (H_i is not 1-dismantlable: exact search),  ? (undecided within the limits);
 *   K1 = the sequence is K-null; >L = |K(H_{i-1})| exceeded the limit L (computation stopped).
 *   With -t a last column gives the CPU time in milliseconds.
 *
 * Summary mode (-S), for large exhaustive runs: only exceptional disks are printed (some iterate whose core is not
 * 1-dismantlable [NO] or undecided [?], or the limit L exceeded [>L]; with -O also the "open" ones, i.e. not K-null
 * within IT steps), and one summary line is written to stderr:
 *   >S total knull_0 ... knull_IT open overL anyNO anyUndecided maxK maxCore
 * where knull_i counts the disks whose first 0-dismantlable iterate is K^i(B).
 *
 * Usage: kdisk [-i IT=5] [-L LIMIT=2000] [-x EXACTMAX=96] [-r GREEDY_TRIES=8] [-s SEED] [-t] [-S [-O]] < disks.g6 > out.tsv
 * Build: gcc -O3 -march=native -o kdisk kdisk.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

typedef uint64_t u64;
#define POP(x) __builtin_popcountll(x)
#define CTZ(x) __builtin_ctzll(x)

typedef struct { int n, W; u64 *a; } Graph;       /* adjacency rows (no loops), W = words per row */

static int DUMP = 0; static int IT = 5, LIMIT = 2000, EXACTMAX = 96, GTRIES = 8, TIMING = 0, SUMMARY = 0, PRINTOPEN = 0;
static u64 rng_state = 88172645463325252ULL;
static u64 rnd(void) { rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17; return rng_state; }

static Graph gnew(int n) { Graph g; g.n = n; g.W = (n + 63) / 64; if (g.W == 0) g.W = 1; g.a = calloc((size_t)n * g.W + 1, sizeof(u64)); return g; }
static inline u64 *row(const Graph *g, int i) { return g->a + (size_t)i * g->W; }
static inline int has(const u64 *s, int j) { return (s[j >> 6] >> (j & 63)) & 1; }
static inline void setb(u64 *s, int j) { s[j >> 6] |= 1ULL << (j & 63); }
static inline void clrb(u64 *s, int j) { s[j >> 6] &= ~(1ULL << (j & 63)); }
static void addedge(Graph *g, int i, int j) { setb(row(g, i), j); setb(row(g, j), i); }
static int popc(const u64 *s, int W) { int c = 0; for (int k = 0; k < W; k++) c += POP(s[k]); return c; }
static long nedges(const Graph *g) { long e = 0; for (int i = 0; i < g->n; i++) e += popc(row(g, i), g->W); return e / 2; }

/* ---------- graph6 ---------- */
static int parse_g6(const char *s, Graph *g) {
    int n; const unsigned char *p = (const unsigned char *)s;
    if (*p == '>') p += 10;                                   /* >>graph6<< header */
    if (*p == 126) { n = ((p[1] - 63) << 12) | ((p[2] - 63) << 6) | (p[3] - 63); p += 4; } else n = *p++ - 63;
    *g = gnew(n);
    int k = 6, x = 0;
    for (int j = 1; j < n; j++) for (int i = 0; i < j; i++) {
        if (k == 6) { x = *p++ - 63; k = 0; }
        if ((x >> (5 - k)) & 1) addedge(g, i, j);
        k++;
    }
    return n;
}


static void write_g6(const Graph *g, FILE *f) {
    int n = g->n;
    if (n < 63) fputc(n + 63, f); else { fputc(126, f); fputc(((n >> 12) & 63) + 63, f); fputc(((n >> 6) & 63) + 63, f); fputc((n & 63) + 63, f); }
    int k = 0, x = 0;
    for (int j = 1; j < n; j++) for (int i = 0; i < j; i++) {
        x = (x << 1) | has(row(g, j), i); if (++k == 6) { fputc(x + 63, f); k = 0; x = 0; }
    }
    if (k) { x <<= (6 - k); fputc(x + 63, f); }
}

/* ---------- 0-core: delete dominated vertices (queue; only neighbors of a deleted vertex need rechecking) ---------- */
static int dominated(const Graph *g, const u64 *alive, int x, u64 *tmp) {
    int W = g->W; const u64 *rx = row(g, x);
    for (int k = 0; k < W; k++) {
        u64 m = rx[k] & alive[k];
        while (m) {
            int y = (k << 6) + CTZ(m); m &= m - 1;
            const u64 *ry = row(g, y); int ok = 1;                 /* N[x] cap alive  subset of  N[y] */
            for (int t = 0; t < W && ok; t++) {
                u64 nx = (rx[t] & alive[t]) | (((x >> 6) == t) ? (1ULL << (x & 63)) : 0);
                u64 ny = ry[t] | (((y >> 6) == t) ? (1ULL << (y & 63)) : 0);
                if (nx & ~ny) ok = 0;
            }
            if (ok) return 1;
        }
    }
    (void)tmp; return 0;
}
/* computes the 0-core of g restricted to alive (alive is modified) and returns the number of remaining vertices */
static int core_mask(const Graph *g, u64 *alive) {
    int n = g->n, W = g->W, cnt = popc(alive, W);
    int *q = malloc(sizeof(int) * (n + 1) * 4); char *inq = calloc(n, 1); int qh = 0, qt = 0, qcap = (n + 1) * 4;
    for (int i = 0; i < n; i++) if (has(alive, i)) { q[qt++] = i; inq[i] = 1; }
    while (qh != qt && cnt > 1) {
        int x = q[qh]; qh = (qh + 1) % qcap; inq[x] = 0;
        if (!has(alive, x)) continue;
        if (dominated(g, alive, x, NULL)) {
            clrb(alive, x); cnt--;
            const u64 *rx = row(g, x);
            for (int k = 0; k < W; k++) { u64 m = rx[k] & alive[k];
                while (m) { int y = (k << 6) + CTZ(m); m &= m - 1; if (!inq[y]) { q[qt] = y; qt = (qt + 1) % qcap; inq[y] = 1; } } }
        }
    }
    free(q); free(inq); return cnt;
}
static Graph induced(const Graph *g, const u64 *mask) {
    int n = g->n, *id = malloc(sizeof(int) * n), m = 0;
    for (int i = 0; i < n; i++) id[i] = has(mask, i) ? m++ : -1;
    Graph h = gnew(m);
    for (int i = 0; i < n; i++) if (id[i] >= 0) {
        const u64 *r = row(g, i);
        for (int k = 0; k < g->W; k++) { u64 x = r[k] & mask[k];
            while (x) { int j = (k << 6) + CTZ(x); x &= x - 1; if (j > i) addedge(&h, id[i], id[j]); } }
    }
    free(id); return h;
}
static Graph core0(const Graph *g) {
    u64 *alive = calloc(g->W, sizeof(u64)); for (int i = 0; i < g->n; i++) setb(alive, i);
    core_mask(g, alive); Graph h = induced(g, alive); free(alive); return h;
}

/* ---------- cliques (Bron--Kerbosch with Tomita pivoting) and clique graph ---------- */
typedef struct { u64 *c; int m, cap, W, overflow, limit; } Cliques;
static void bk(const Graph *g, u64 *R, u64 *P, u64 *X, Cliques *C) {
    int W = g->W, emptyP = 1, emptyX = 1;
    for (int k = 0; k < W; k++) { if (P[k]) emptyP = 0; if (X[k]) emptyX = 0; }
    if (emptyP) {
        if (emptyX) {
            if (C->m >= C->limit) { C->overflow = 1; return; }
            if (C->m == C->cap) { C->cap *= 2; C->c = realloc(C->c, sizeof(u64) * (size_t)C->cap * W); }
            memcpy(C->c + (size_t)C->m * W, R, sizeof(u64) * W); C->m++;
        }
        return;
    }
    int best = -1, bc = -1;                                   /* pivot u in P cup X maximizing |P cap N(u)| */
    for (int k = 0; k < W; k++) { u64 m = P[k] | X[k];
        while (m) { int u = (k << 6) + CTZ(m); m &= m - 1; int c = 0; const u64 *ru = row(g, u);
            for (int t = 0; t < W; t++) c += POP(P[t] & ru[t]); if (c > bc) { bc = c; best = u; } } }
    u64 *cand = malloc(sizeof(u64) * W), *R2 = malloc(sizeof(u64) * W), *P2 = malloc(sizeof(u64) * W), *X2 = malloc(sizeof(u64) * W);
    const u64 *rp = row(g, best);
    for (int k = 0; k < W; k++) cand[k] = P[k] & ~rp[k];
    for (int k = 0; k < W && !C->overflow; k++) {
        while (cand[k] && !C->overflow) {
            int v = (k << 6) + CTZ(cand[k]); cand[k] &= cand[k] - 1;
            const u64 *rv = row(g, v);
            for (int t = 0; t < W; t++) { R2[t] = R[t]; P2[t] = P[t] & rv[t]; X2[t] = X[t] & rv[t]; }
            setb(R2, v);
            bk(g, R2, P2, X2, C);
            clrb(P, v); setb(X, v);
        }
    }
    free(cand); free(R2); free(P2); free(X2);
}
/* returns the clique graph, or a graph with n = -1 if there are more than limit cliques */
static Graph clique_graph(const Graph *g, int limit) {
    int W = g->W;
    Cliques C = { malloc(sizeof(u64) * 64 * W), 0, 64, W, 0, limit };
    u64 *R = calloc(W, 8), *P = calloc(W, 8), *X = calloc(W, 8);
    for (int i = 0; i < g->n; i++) setb(P, i);
    bk(g, R, P, X, &C);
    free(R); free(P); free(X);
    Graph k;
    if (C.overflow) { k.n = -1; k.W = 0; k.a = NULL; free(C.c); return k; }
    k = gnew(C.m);
    for (int p = 0; p < C.m; p++) for (int q = p + 1; q < C.m; q++) {
        const u64 *a = C.c + (size_t)p * W, *b = C.c + (size_t)q * W;
        for (int t = 0; t < W; t++) if (a[t] & b[t]) { addedge(&k, p, q); break; }
    }
    free(C.c); return k;
}

/* ---------- 1-dismantlability ---------- */
/* v is 1-dominated in g[alive]  <=>  N(v) cap alive is nonempty and 0-dismantlable */
static int one_dominated(const Graph *g, const u64 *alive, int v) {
    int W = g->W; const u64 *rv = row(g, v);
    u64 *mask = malloc(sizeof(u64) * W); int d = 0;
    for (int k = 0; k < W; k++) { mask[k] = rv[k] & alive[k]; d += POP(mask[k]); }
    int res;
    if (d == 0) res = 0;
    else if (d <= 2) res = 1;                                  /* K1, K2 or 2K1? 2K1 is not 0-dismantlable */
    else res = -1;
    if (d == 2) {                                             /* two vertices: 0-dismantlable iff adjacent */
        int u[2], c = 0;
        for (int k = 0; k < W; k++) { u64 m = mask[k]; while (m) { u[c++] = (k << 6) + CTZ(m); m &= m - 1; } }
        res = has(row(g, u[0]), u[1]);
    }
    if (res < 0) { Graph h = induced(g, mask); u64 *al = calloc(h.W, 8); for (int i = 0; i < h.n; i++) setb(al, i);
        res = core_mask(&h, al) == 1; free(al); free(h.a); }
    free(mask); return res;
}
/* randomized greedy deletion of 1-dominated vertices; statuses are cached and recomputed only for the
   neighbors of a deleted vertex (the neighborhoods of the other vertices do not change) */
static int greedy1(const Graph *g, int tries) {
    int n = g->n, W = g->W;
    if (n <= 1) return 1;
    u64 *alive = malloc(sizeof(u64) * W); signed char *st = malloc(n); int *ord = malloc(sizeof(int) * n);
    int success = 0;
    for (int t = 0; t < tries && !success; t++) {
        memset(alive, 0, sizeof(u64) * W); for (int i = 0; i < n; i++) setb(alive, i);
        for (int i = 0; i < n; i++) st[i] = (signed char)one_dominated(g, alive, i);
        int cnt = n;
        while (cnt > 1) {
            int start = (int)(rnd() % n), v = -1;
            for (int s = 0; s < n; s++) { int x = (start + s) % n; if (has(alive, x) && st[x]) { v = x; break; } }
            if (v < 0) break;
            clrb(alive, v); cnt--;
            const u64 *rv = row(g, v);
            for (int k = 0; k < W; k++) { u64 m = rv[k] & alive[k];
                while (m) { int u = (k << 6) + CTZ(m); m &= m - 1; st[u] = (signed char)one_dominated(g, alive, u); } }
        }
        if (cnt == 1) success = 1;
        (void)ord;
    }
    free(alive); free(st); free(ord); return success;
}
/* exact search: depth-first over the subsets reachable by deleting 1-dominated vertices, with a hash set of
   subsets already known to fail; returns 1, 0, or -1 (node limit reached) */
typedef struct { u64 *keys; size_t sz, cnt; int W; } HSet;
static u64 hkey(const u64 *s, int W) { u64 h = 0x9E3779B97F4A7C15ULL; for (int k = 0; k < W; k++) { h ^= s[k]; h *= 0xff51afd7ed558ccdULL; h ^= h >> 33; } return h; }
static int hs_has(HSet *H, const u64 *s) { size_t p = hkey(s, H->W) & (H->sz - 1);
    while (1) { u64 *k = H->keys + p * (H->W + 1); if (!k[0]) return 0; if (!memcmp(k + 1, s, 8 * H->W)) return 1; p = (p + 1) & (H->sz - 1); } }
static void hs_add(HSet *H, const u64 *s);
static void hs_grow(HSet *H) { HSet N = { calloc(H->sz * 2 * (H->W + 1), 8), H->sz * 2, 0, H->W };
    for (size_t i = 0; i < H->sz; i++) { u64 *k = H->keys + i * (H->W + 1); if (k[0]) hs_add(&N, k + 1); }
    free(H->keys); *H = N; }
static void hs_add(HSet *H, const u64 *s) { if (2 * (H->cnt + 1) > H->sz) hs_grow(H);
    size_t p = hkey(s, H->W) & (H->sz - 1);
    while (1) { u64 *k = H->keys + p * (H->W + 1); if (!k[0]) { k[0] = 1; memcpy(k + 1, s, 8 * H->W); H->cnt++; return; } p = (p + 1) & (H->sz - 1); } }
static long nodes, NODELIMIT = 2000000;
static int ex_rec(const Graph *g, u64 *S, int cnt, HSet *fail) {
    if (cnt <= 1) return 1;
    if (hs_has(fail, S)) return 0;
    if (++nodes > NODELIMIT) return -1;
    int W = g->W;
    u64 *T = malloc(8 * W);
    for (int k = 0; k < W; k++) { u64 m = S[k];
        while (m) { int v = (k << 6) + CTZ(m); m &= m - 1;
            if (one_dominated(g, S, v)) {
                memcpy(T, S, 8 * W); clrb(T, v);
                int r = ex_rec(g, T, cnt - 1, fail);
                if (r != 0) { free(T); return r; }
            } } }
    free(T); hs_add(fail, S); return 0;
}
static int exact1(const Graph *g) {
    if (g->n > EXACTMAX) return -1;
    HSet fail = { calloc(1024 * (g->W + 1), 8), 1024, 0, g->W };
    u64 *S = calloc(g->W, 8); for (int i = 0; i < g->n; i++) setb(S, i);
    nodes = 0; int r = ex_rec(g, S, g->n, &fail);
    free(S); free(fail.keys); return r;
}
static const char *status1(const Graph *h) {
    if (h->n == 1) return "0";
    if (greedy1(h, GTRIES)) return "1";
    int r = exact1(h);
    return r == 1 ? "1" : (r == 0 ? "NO" : "?");
}

int main(int argc, char **argv) {
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-i")) IT = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-L")) LIMIT = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-x")) EXACTMAX = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-r")) GTRIES = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-s")) rng_state = strtoull(argv[++i], NULL, 10) | 1;
        else if (!strcmp(argv[i], "-t")) TIMING = 1;
        else if (!strcmp(argv[i], "-S")) SUMMARY = 1;
        else if (!strcmp(argv[i], "-O")) PRINTOPEN = 1;
        else if (!strcmp(argv[i], "-D")) DUMP = 1;
        else { fprintf(stderr, "unknown option %s\n", argv[i]); return 1; }
    }
    char line[1 << 16];
    static char out[1 << 16];
    long long total = 0, knull[64] = {0}, nopen = 0, nbig = 0, nno = 0, nunk = 0; int maxK = 0, maxC = 0;
    while (fgets(line, sizeof line, stdin)) {
        size_t L = strlen(line); while (L && (line[L - 1] == '\n' || line[L - 1] == '\r')) line[--L] = 0;
        if (!L) continue;
        clock_t t0 = clock();
        int o = 0, big = 0, hasno = 0, hasunk = 0, first0 = -1;
        Graph B; parse_g6(line, &B);
        o += sprintf(out + o, "%s\t%d\t%ld", line, B.n, nedges(&B));
        Graph H = core0(&B);
        const char *st = status1(&H);
        if (!strcmp(st, "NO")) hasno = 1; if (!strcmp(st, "?")) hasunk = 1;
        if (H.n > maxC) maxC = H.n;
        o += sprintf(out + o, "\t%d:%d:%s", B.n, H.n, st);
        free(B.a);
        if (H.n == 1) first0 = 0;
        for (int i = 1; i <= IT && first0 < 0; i++) {
            Graph K = clique_graph(&H, LIMIT);
            if (K.n < 0) { o += sprintf(out + o, "\t>%d", LIMIT); big = 1; if (DUMP) { fprintf(stderr, "CORE%d\t%s\t%d\t", i - 1, line, H.n); write_g6(&H, stderr); fputc(10, stderr); } break; }
            if (K.n > maxK) maxK = K.n;
            Graph H2 = core0(&K);
            st = status1(&H2);
            if (!strcmp(st, "NO")) hasno = 1; if (!strcmp(st, "?")) hasunk = 1;
            if (H2.n > maxC) maxC = H2.n;
            o += sprintf(out + o, "\t%d:%d:%s", K.n, H2.n, st);
            free(K.a); free(H.a); H = H2;
            if (H.n == 1) first0 = i;
        }
        if (first0 >= 0) o += sprintf(out + o, "\tK1");
        free(H.a);
        if (TIMING) o += sprintf(out + o, "\t%.1f", 1000.0 * (clock() - t0) / CLOCKS_PER_SEC);
        total++;
        if (first0 >= 0) knull[first0]++; else if (big) nbig++; else nopen++;
        nno += hasno; nunk += hasunk;
        int show = !SUMMARY || hasno || hasunk || big || (PRINTOPEN && first0 < 0);
        if (show) { fputs(out, stdout); fputc('\n', stdout); fflush(stdout); }
    }
    if (SUMMARY) {
        fprintf(stderr, ">S %lld", total);
        for (int i = 0; i <= IT; i++) fprintf(stderr, " %lld", knull[i]);
        fprintf(stderr, " %lld %lld %lld %lld %d %d\n", nopen, nbig, nno, nunk, maxK, maxC);
    }
    return 0;
}
