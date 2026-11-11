/* gen_disks.c -- exhaustive generation of Whitney triangulations of the disk (flag disks) up to isomorphism.
 *
 * Every flag disk D with n >= 4 vertices has a boundary vertex v whose neighborhood is a path x_0...x_k of
 * boundary vertices of D - v, with D - v again a flag disk (Proposition 2.7 of the paper).  Conversely, if D' is a
 * flag disk and x_0...x_k (k >= 1) is a path of consecutive boundary vertices of D' that is an induced path of D',
 * then adding a vertex v adjacent to x_0,...,x_k gives a flag disk.  Hence all flag disks on n vertices are
 * obtained from those on n-1 vertices by this augmentation; isomorphic copies are removed with nauty's canonical
 * labelling and a hash table.
 *
 * Usage:  gen_disks NMAX [outdir]
 *   writes outdir/disks_<n>.g6 (graph6, one disk per line) for 3 <= n <= NMAX and prints the counts and times.
 * Build:  gcc -O3 -march=native -o gen_disks gen_disks.c -DWORDSIZE=64 -DMAXN=WORDSIZE -lnautyL1
 *         (or link the nauty sources compiled with MAXN=WORDSIZE).  NMAX <= 63.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "nauty.h"
#include "gtools.h"

typedef struct { setword *rows; size_t cnt, cap; int n; size_t *tab; size_t tabsz; } Level;

static unsigned long long hashrows(const setword *r, int n) {
    unsigned long long h = 1469598103934665603ULL;
    for (int i = 0; i < n; i++) { h ^= (unsigned long long)r[i]; h *= 1099511628211ULL; h ^= h >> 29; }
    return h;
}
static void level_init(Level *L, int n) {
    L->n = n; L->cnt = 0; L->cap = 1024; L->rows = malloc(L->cap * n * sizeof(setword));
    L->tabsz = 4096; L->tab = calloc(L->tabsz, sizeof(size_t));          /* stores index+1, 0 = empty */
}
static void level_rehash(Level *L) {
    free(L->tab); L->tabsz *= 2; L->tab = calloc(L->tabsz, sizeof(size_t));
    for (size_t i = 0; i < L->cnt; i++) {
        size_t p = hashrows(L->rows + i * L->n, L->n) & (L->tabsz - 1);
        while (L->tab[p]) p = (p + 1) & (L->tabsz - 1);
        L->tab[p] = i + 1;
    }
}
static int level_insert(Level *L, const setword *r) {           /* returns 1 if new */
    int n = L->n;
    size_t p = hashrows(r, n) & (L->tabsz - 1);
    while (L->tab[p]) {
        if (memcmp(L->rows + (L->tab[p] - 1) * n, r, n * sizeof(setword)) == 0) return 0;
        p = (p + 1) & (L->tabsz - 1);
    }
    if (L->cnt == L->cap) { L->cap *= 2; L->rows = realloc(L->rows, L->cap * n * sizeof(setword)); }
    memcpy(L->rows + L->cnt * n, r, n * sizeof(setword));
    L->tab[p] = ++L->cnt;
    if (2 * L->cnt > L->tabsz) level_rehash(L);
    return 1;
}

static void canon_insert(Level *L, graph *g, int n) {
    int lab[MAXN], ptn[MAXN], orbits[MAXN];
    graph cg[MAXN];
    DEFAULTOPTIONS_GRAPH(options); statsblk stats;
    options.getcanon = TRUE;
    densenauty(g, lab, ptn, orbits, &options, &stats, 1, n, cg);
    level_insert(L, cg);
}

/* boundary cycle of a flag disk given by rows g (n vertices); returns its length b, cycle in cyc[] */
static int boundary_cycle(const setword *g, int n, int *cyc) {
    int bdeg[MAXN]; setword bnb[MAXN];
    for (int i = 0; i < n; i++) { bnb[i] = 0; bdeg[i] = 0; }
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (ISELEMENT(&g[i], j) && POPCOUNT(g[i] & g[j]) == 1) {   /* edge in exactly one triangle */
                bnb[i] |= bit[j]; bnb[j] |= bit[i]; bdeg[i]++; bdeg[j]++;
            }
    int start = -1;
    for (int i = 0; i < n; i++) if (bdeg[i]) { start = i; break; }
    if (start < 0) return 0;
    int b = 0, prev = -1, cur = start;
    do {
        cyc[b++] = cur;
        int nxt = -1;
        for (int j = 0; j < n; j++) if (ISELEMENT(&bnb[cur], j) && j != prev) { nxt = j; break; }
        prev = cur; cur = nxt;
    } while (cur != start && cur >= 0 && b <= n);
    return b;
}

static void augment(const setword *g, int n, Level *next) {
    int cyc[MAXN]; int b = boundary_cycle(g, n, cyc);
    graph h[MAXN];
    for (int s = 0; s < b; s++)
        for (int k = 1; k <= b - 1; k++) {
            int ok = 1;
            for (int i = 0; i <= k && ok; i++)
                for (int j = i + 2; j <= k; j++)
                    if (ISELEMENT(&g[cyc[(s + i) % b]], cyc[(s + j) % b])) { ok = 0; break; }
            if (!ok) break;                                   /* longer paths starting at s contain this chord */
            for (int i = 0; i < n; i++) h[i] = g[i];
            h[n] = 0;
            for (int i = 0; i <= k; i++) { int x = cyc[(s + i) % b]; h[x] |= bit[n]; h[n] |= bit[x]; }
            canon_insert(next, h, n + 1);
        }
}

int main(int argc, char **argv) {
    int NMAX = argc > 1 ? atoi(argv[1]) : 10;
    const char *dir = argc > 2 ? argv[2] : ".";
    if (NMAX > MAXN - 1) { fprintf(stderr, "NMAX too large\n"); return 1; }
    nauty_check(WORDSIZE, 1, MAXN, NAUTYVERSIONID);
    Level cur; level_init(&cur, 3);
    graph t[3] = {bit[1] | bit[2], bit[0] | bit[2], bit[0] | bit[1]};
    canon_insert(&cur, t, 3);
    for (int n = 3; n <= NMAX; n++) {
        clock_t t0 = clock();
        char fn[512]; snprintf(fn, sizeof fn, "%s/disks_%d.g6", dir, n);
        FILE *f = fopen(fn, "w");
        for (size_t i = 0; i < cur.cnt; i++) fputs(ntog6(cur.rows + i * n, 1, n), f);
        fclose(f);
        Level nxt; level_init(&nxt, n + 1);
        if (n < NMAX) for (size_t i = 0; i < cur.cnt; i++) augment(cur.rows + i * n, n, &nxt);
        printf("n=%2d  disks=%zu  (augmentation to n+1: %.2f s)\n", n, cur.cnt, (double)(clock() - t0) / CLOCKS_PER_SEC);
        fflush(stdout);
        free(cur.rows); free(cur.tab); cur = nxt;
    }
    return 0;
}
