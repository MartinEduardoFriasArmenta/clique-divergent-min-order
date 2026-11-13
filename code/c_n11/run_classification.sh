#!/bin/bash
# Theorem 4.4: K-behavior of the 1-dismantlable, non-0-dismantlable graphs on 11 vertices.
# Input: the 2000 parts written by correr_c.sh (lines "graph6<TAB>|K(G)|<TAB>state", state 0, 1 or NO)
# in resultados_n11/.  Output: stage files like those in ../../results/classification_n11/.
#   state 0 : K(G) 0-dismantlable           -> K-null in one step                  12,240,034 graphs
#   state NO: K(G) not 1-dismantlable       -> the 81 K-divergent graphs of Theorem 4.3
#   state 1 : K(G) 1-dismantlable, not 0-d. -> iterate K with 0-cores (kdisk)       3,653,602 graphs
set -e
cd "$(dirname "$0")"
gcc -O3 -march=native -o kdisk ../exhaustive_disks_C/kdisk.c
# graphs with state 1 -> state1.g6 (the parts are in resultados_n11/, written by correr_c.sh)
for f in resultados_n11/parte_*.tsv.gz; do zcat "$f"; done | awk -F'\t' '$3=="1"{print $1}' > state1.g6
split -n l/2 -d state1.g6 half_                      # two halves, one per core
for h in 00 01; do ./kdisk -i 12 -L 3000 -S -O < half_$h > stage1_$h.out 2> stage1_$h.sum & done; wait
cat stage1_00.out stage1_01.out | cut -f1 > over.g6  # graphs whose iterates exceeded |K(H)| = 3000 (18 graphs)
./kdisk -i 12 -L 40000 -t < over.g6 > stage2_L40000.tsv      # 8 of them are K-null; 10 remain
awk -F'\t' '$(NF-1)==">40000"{print $1}' stage2_L40000.tsv > ten.g6
./kdisk -D -i 12 -L 150000 -t < ten.g6 > stage3_L150000.tsv 2> cores.txt   # 1 more is K-null (K^6); 9 remain (memory 150000^2/8 bytes)
