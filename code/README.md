# Programs for "A 1-dismantlable, collapsible, clique-divergent graph of minimum order"

Requirements: Python 3 (tested with 3.13), `networkx`, `numpy`; for the exhaustive enumerations, `geng` from nauty 2.8.x; a C compiler for `c_n11/`.
(Some comments and console messages inside the older programs are in Spanish; `README_scripts.md` is the original Spanish version of this file.)

## Auxiliary modules
- `discoid.py` — bitmask routines: `zero_dominated`, `zero_dismantlable`, `one_dominated`, class `Discoid` (exact decision of 1-dismantlability by memoized search over vertex subsets), `is_discoid(G, with_order=True)` (returns a deletion order), `clique_graph(G)`.
- `exp15.py` — graph6 reader `parse_g6` and Bron–Kerbosch clique enumeration `cliques`.
- `homol.py` — `betti_z2(G, maxdim)`: reduced Betti numbers over Z_2 of the clique complex; `python3 homol.py <graph6>` prints them.
- `retract.py` — induced octahedra `octahedra(G, k)` and backtracking search for a retraction `retraction(G, S)`.

## Experiments (section numbers refer to the paper)
1. `propiedades_A.py [graph6]` — Theorems 3.1–3.3 and Proposition 5.1 (default: A = `J~Z}eY~qzd_`; with `JrqipkN~~~?`, the graph of Larrión–Pizaña–Villarroel-Flores). ~15 s / ~45 s.
2. `exh_disc2.py <geng> <n> <out.tsv>` — Theorem 4.1 (n = 10 takes 1–2 h on one core). Outputs in `../data/discoides_{8,9,10}.tsv`.
3. `kbehav.py <discoides_n.tsv>` — Theorem 4.2 (K-nullity by iterating K and taking 0-cores).
4. `c_n11/` — Theorem 4.3: `k1proc.c` (OUTPROC module for geng), `compilar.sh`, `correr_c.sh` (2000 resumable parts; ~5 core-hours), `analisis_final.py`. Result: `../data/contraejemplos_n11.tsv` (the 81 graphs).
4b. `c_n11/run_classification.sh` — Theorem 4.4: K-behavior of the other 1-dismantlable graphs on 11 vertices, with `exhaustive_disks_C/kdisk.c` (about 10 minutes for the first stage, plus the larger limits). Results: `../results/classification_n11/`.
5. `certificado_K2/` — Theorem 4.3, the graph G_* = `JEhfrz{~~~?`: `retract_oct.py` (retraction onto an octahedron as a constraint-satisfaction problem), `k2_oct2.py` (search for the O_16 retract of K^2(G_*)), `verificar_k2.py` (independent verification). Certificate: `../data/k2_retraccion_JEhfrz.json`.
6. `susp11.py <geng> <out.tsv>` — Proposition 5.2 (the 37 suspensions on 11 vertices). Output `../data/susp11.tsv` (~7 min).
7. `family_Gl/` — Section 6, the family G_l (standalone folder):
   - `python3 family_Gl.py [L]` checks, for l = 1..L (default 16; a few seconds): the deletion order of Theorem 6.3(1); flag sphere, no 0-dominated vertices, not a suspension, no clique separator; the explicit retraction r_l of G_l onto Susp(D_l) ≅ G_LPV (Theorem 6.3(3)); and, independently, retractions of K(G_l) onto an induced O_4 (CSP certificates for l ≤ 8; uniform retraction plus locality check for l ≥ 9). Certificates: `../data/family_Gl_certificates.json`.
   - `python3 draw_Gl.py` draws Figure 4 (`fig_Gl.pdf`).
8. Section 7 (exhaustive, n ≤ 18): `exhaustive_disks_C/` (see its README); results in `../results/exhaustive_disks_n18/`. Preliminary Python experiment: `disk_exp.py <geng> <minutes> <out.tsv> [seed]`, data `../data/disk_results.tsv` (1403 disks).

Quick reproduction of the main example: `python3 propiedades_A.py`.
