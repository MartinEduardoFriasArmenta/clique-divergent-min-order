# A 1-dismantlable, collapsible, clique-divergent graph of minimum order

Code, data and computational certificates accompanying the paper

> J.F. Espinoza-Fierro, M.E. Frías-Armenta, H.A. Hernández Hernández,
> *A 1-dismantlable, collapsible, clique-divergent graph of minimum order*, 2026 (submitted).

The paper presents the graph `A` (graph6 `J~Z}eY~qzd_`), an 11-vertex graph that is 1-dismantlable (hence
nonevasive and collapsible) but clique-divergent; proves that 11 is the minimum order of such a graph;
classifies the 81 examples of minimum order; constructs an infinite family on flag 2-spheres; and verifies the
disk conjecture of Larrión, Neumann-Lara and Pizaña for all 550,551,244 Whitney triangulations of the disk with at
most 18 vertices.

## Layout

| Directory | Contents |
|---|---|
| `code/` | All programs (Python 3 + C). See `code/README.md` for a description of each one and the theorem it supports. |
| `data/` | Outputs of the enumerations: the 81 graphs on 11 vertices (`contraejemplos_n11.tsv`), the 37 suspensions, the 1-dismantlable graphs on 8–10 vertices, the preliminary disk experiment, the retraction certificate for `K^2(G_*)`, the certificates for the family `G_l`. |
| `results/exhaustive_disks_n18/` | Summary of the exhaustive run over all Whitney disks with n ≤ 18 vertices (Section 7): per-n tables in TSV and LaTeX, the 39 hard disks, raw part summaries. |

## Quick reproduction

```bash
pip install networkx numpy            # Python 3.10+
sudo apt install nauty libnauty-dev   # geng, labelg, and nauty.h (or: brew install nauty)

cd code
python3 propiedades_A.py              # Theorems 3.1–3.3: A, K(A), the O_4 retraction (~15 s)
python3 propiedades_A.py 'JrqipkN~~~?' # Proposition 5.1: the graph of Larrión–Pizaña–Villarroel-Flores
python3 exh_disc2.py geng 9 out9.tsv  # Theorem 4.1 for n = 9 (n = 10 takes 1–2 h)
cd family_Gl && python3 family_Gl.py  # Section 6: the family G_l, l = 1..16 (seconds)
cd ../exhaustive_disks_C && make && ./run_exhaustive.sh 3 14 2 && python3 summarize.py out   # Section 7 up to n = 14 (minutes)
```

The full computations of the paper take about 5 core-hours (Theorem 4.3, `code/c_n11/`) and about 7 hours on
two cores (Theorem 7.2, `code/exhaustive_disks_C/`, up to n = 18).

## Theorem → program

| Result | Program(s) | Data |
|---|---|---|
| Thm 3.1–3.3, Prop 5.1 | `propiedades_A.py` | `data/colapsos_*.txt` |
| Thm 4.1, 4.2 | `exh_disc2.py`, `kbehav.py` | `data/discoides_{8,9,10}.tsv` |
| Thm 4.3 | `c_n11/` (`k1proc.c` + geng), `analisis_final.py` | `data/contraejemplos_n11.tsv` |
| Thm 4.3, graph `G_*` | `certificado_K2/` | `data/k2_retraccion_JEhfrz.json` |
| Prop 5.2 | `susp11.py` | `data/susp11.tsv` |
| Thm 6.3 | `family_Gl/family_Gl.py` | `data/family_Gl_certificates.json` |
| Thm 7.2 | `exhaustive_disks_C/` | `results/exhaustive_disks_n18/` |

## Citation

See `CITATION.cff`. Code is released under the MIT License; data and certificates under CC BY 4.0.
