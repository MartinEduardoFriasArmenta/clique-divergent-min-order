/* Plug-in para geng (nauty): para cada gráfica conexa G que produce geng decide
 *   - si G es 0-desmantelable (borrado de vértices dominados hasta K1),
 *   - si G es 1-desmantelable (búsqueda exacta con memo sobre subconjuntos),
 *   - para las 1-desmantelables no 0-desmantelables: calcula K(G) y decide si K(G) es
 *     0-desmantelable, 1-desmantelable, NO 1-desmantelable, o indecidido (tope de estados).
 * Escribe en la salida de geng una línea "graph6<TAB>|K(G)|<TAB>estado" por cada G
 * 1-desmantelable no 0-desmantelable (estado: 0, 1, NO, ?), y al final, en stderr, los conteos
 * con el mismo formato que exh_disc2.py.
 * Compilación (ver compilar.sh):
 *   gcc -O3 -march=native -DMAXN=WORDSIZE -DWORDSIZE=32 -DOUTPROC=k1out -DSUMMARY=k1sum \
 *       -o geng_k1 geng.c k1proc.c gtools.c nauty.c nautil.c naugraph.c schreier.c naurng.c
 */
#include "gtools.h"
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

typedef uint64_t u64;
static int lastn = 0;
static long long c_tot = 0, c_disc = 0, c_kndz = 0, c_bad = 0, c_und = 0;
#define POP(x) __builtin_popcountll(x)
#define LOW(x) __builtin_ctzll(x)

/* ---------- 0-desmantelabilidad (máscaras de 64 bits) ---------- */
static inline int zero_dominated(const u64 *adj, u64 S, int v) {
    u64 Nv = (adj[v] | (1ULL << v)) & S, m = adj[v] & S;
    while (m) { int w = LOW(m); m &= m - 1;
        if ((Nv & ~((adj[w] | (1ULL << w)) & S)) == 0) return 1; }
    return 0;
}
static int zero_dismantlable(const u64 *adj, u64 S) {
    if (!S) return 0;
    int changed = 1;
    while (changed && POP(S) > 1) {
        changed = 0; u64 m = S;
        while (m) { int v = LOW(m); m &= m - 1;
            if (zero_dominated(adj, S, v)) { S &= ~(1ULL << v); changed = 1; break; } }
    }
    return POP(S) == 1;
}
static inline int one_dominated(const u64 *adj, u64 S, int v) {
    u64 N = adj[v] & S; return N && zero_dismantlable(adj, N);
}

/* ---------- 1-desmantelabilidad exacta: DFS con memo de fracasos (tabla hash) ---------- */
#define HBITS 23
#define HSIZE (1u << HBITS)
static u64 *htab = NULL;          /* claves = conjunto+1 (0 = vacío) */
static u64 *hused = NULL; static long nused = 0;   /* para limpiar rápido */
static long cap = 3000000;        /* tope de estados */
static int aborted;
static int hfind(u64 key) {
    u64 k = key + 1; unsigned h = (unsigned)((k * 0x9E3779B97F4A7C15ULL) >> (64 - HBITS));
    while (htab[h]) { if (htab[h] == k) return 1; h = (h + 1) & (HSIZE - 1); }
    return 0;
}
static void hadd(u64 key) {
    u64 k = key + 1; unsigned h = (unsigned)((k * 0x9E3779B97F4A7C15ULL) >> (64 - HBITS));
    while (htab[h]) { if (htab[h] == k) return; h = (h + 1) & (HSIZE - 1); }
    htab[h] = k; hused[nused++] = h;
}
static void hclear(void) { for (long i = 0; i < nused; i++) htab[hused[i]] = 0; nused = 0; }

static int dfs(const u64 *adj, u64 S) {
    if (POP(S) == 1) return 1;
    if (hfind(S)) return 0;
    if (nused >= cap) { aborted = 1; return 0; }
    int c0[64], c1[64], n0 = 0, n1 = 0; u64 m = S;
    while (m) { int v = LOW(m); m &= m - 1;
        if (zero_dominated(adj, S, v)) c0[n0++] = v;
        else if (one_dominated(adj, S, v)) c1[n1++] = v; }
    for (int i = 0; i < n0; i++) { if (dfs(adj, S & ~(1ULL << c0[i]))) return 1; if (aborted) return 0; }
    for (int i = 0; i < n1; i++) { if (dfs(adj, S & ~(1ULL << c1[i]))) return 1; if (aborted) return 0; }
    hadd(S);
    return 0;
}
/* 1: sí, 0: no, -1: indecidido */
static int one_dismantlable(const u64 *adj, int n) {
    u64 S = (n == 64) ? ~0ULL : ((1ULL << n) - 1);
    aborted = 0; int r = dfs(adj, S); hclear();
    return aborted ? -1 : r;
}

/* ---------- clanes (Bron–Kerbosch con pivote) ---------- */
static u64 cl[64]; static int ncl;
static void bk(const u64 *adj, u64 R, u64 P, u64 X) {
    if (!P && !X) { if (ncl < 64) cl[ncl] = R; ncl++; return; }
    u64 PX = P | X, m = PX; int piv = LOW(PX), best = -1;
    while (m) { int u = LOW(m); m &= m - 1; int c = POP(P & adj[u]); if (c > best) { best = c; piv = u; } }
    u64 cand = P & ~adj[piv];
    while (cand) { int v = LOW(cand); cand &= cand - 1;
        bk(adj, R | (1ULL << v), P & adj[v], X & adj[v]); P &= ~(1ULL << v); X |= 1ULL << v; }
}

void k1out(FILE *f, graph *g, int n) {
    int m = SETWORDSNEEDED(n);
    u64 adj[64]; c_tot++; lastn = n;
    for (int i = 0; i < n; i++) { set *gi = GRAPHROW(g, i, m); adj[i] = 0;
        for (int j = 0; j < n; j++) if (ISELEMENT(gi, j)) adj[i] |= 1ULL << j; }
    u64 S = (1ULL << n) - 1;
    if (zero_dismantlable(adj, S)) return;
    if (!htab) { htab = calloc(HSIZE, sizeof(u64)); hused = malloc(sizeof(u64) * (cap + 64)); }
    int r = one_dismantlable(adj, n);
    if (r != 1) return;                 /* para n<=11 nunca es indecidido */
    c_disc++;
    ncl = 0; bk(adj, 0, S, 0);
    const char *st; int nk = ncl;
    if (ncl > 64) { st = "?"; c_und++; c_kndz++; }
    else {
        u64 K[64];
        for (int a = 0; a < ncl; a++) { K[a] = 0; for (int b = 0; b < ncl; b++) if (a != b && (cl[a] & cl[b])) K[a] |= 1ULL << b; }
        u64 SK = (ncl == 64) ? ~0ULL : ((1ULL << ncl) - 1);
        if (zero_dismantlable(K, SK)) st = "0";
        else { c_kndz++; int rk = one_dismantlable(K, ncl);
            if (rk == 1) st = "1"; else if (rk == 0) { st = "NO"; c_bad++; } else { st = "?"; c_und++; } }
    }
    char *s = ntog6(g, m, n); size_t L = strlen(s); if (L && s[L-1] == '\n') s[L-1] = 0;
    fprintf(f, "%s\t%d\t%s\n", s, nk, st);
}

void k1sum(nauty_counter nout, double cpu) {
    fprintf(stderr, "n=%d conexas=%lld 1-desm_no_0-desm=%lld K_no_0-desm=%lld K_no_1-desm=%lld indecidido=%lld tiempo=%.0fs\n",
            lastn, c_tot, c_disc, c_kndz, c_bad, c_und, cpu);
}
