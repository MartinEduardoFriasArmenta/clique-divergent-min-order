# K-behavior of the 1-dismantlable graphs on 11 vertices

Complete classification of the 15,893,717 graphs on 11 vertices that are 1-dismantlable and not 0-dismantlable
(Theorem 4.3), by the behavior of the iterated clique graphs. Programs in C (`code/kdisk.c`, the program of
Section 7 with an added option `-D` that prints the last core reached).

| class | graphs |
|---|---|
| K(G) not 1-dismantlable: K-divergent (Theorem 4.3) | 81 |
| K-null, first 0-dismantlable iterate K^1(G) | 12,240,034 |
| K-null, K^2(G) | 3,371,703 |
| K-null, K^3(G) | 267,210 |
| K-null, K^4(G) | 14,507 |
| K-null, K^5(G) | 167 |
| K-null, K^6(G) | 6 |
| **undecided** (K(G) is 1-dismantlable; some iterate has a core whose clique graph has more than 150,000 vertices) | **9** |
| total | 15,893,717 |

Stages (see `run_classification.sh`): the 3,653,602 graphs whose clique graph is 1-dismantlable but not 0-dismantlable
were iterated with `kdisk -i 12 -L 3000` (3,653,584 decided; 18 exceeded the limit; `stage1_*.sum`), the 18 with `-L 40000`
(8 decided, `stage2_L40000.tsv`) and the remaining 10 with `-L 150000` (1 decided, `stage3_L150000.tsv`). The nine
graphs left are in `undecided_9.tsv`, with the sizes |K(H_{i-1})| : |H_i| of the iterations; for all of them every core
computed is 1-dismantlable. Three of them are suspensions (40-41 edges); the core of 159 vertices of one of them has more than six million
maximal cliques.

Consequence for the paper: every 1-dismantlable graph on 11 vertices is K-null, except the 81 divergent ones and
these nine, whose K-behavior is open; so "the exact list of examples of minimum order" is established only up to them. The counts agree with an independent classification (reviewer's report): 12,240,034 / 3,371,703 /
267,210 / 14,507 / 167 / 5 K-null graphs (the sixth K^6 graph, `JCqjfT{|fb?` in graph6 of this run, was left
undecided there) and ten undecided graphs, the same ten up to isomorphism (checked with `labelg`).
