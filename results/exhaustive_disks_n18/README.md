# Exhaustive disk experiment, 3 <= n <= 18 (run 2026-11-08, 06:00-11:02, 2 cores)
550,551,244 flag disks (Whitney triangulations of the disk), all of them up to isomorphism.
Programs: exhaustive_disks_C (gen_disks_ca | kdisk -S -O -i 5 -L 2000).
- exhaustive_summary.tsv / exhaustive_table.tex: per n, number of disks, K-null with index i (first 0-dismantlable iterate K^i(B), i = 0..5), open (not K-null within 5 steps), >L (an iterate K(H_4) with more than 2000 vertices), NO / ? (an iterate whose core is not 1-dismantlable / undecided), max |K(H)|, max |core|.
- limit_exceeded_n18.tsv: the 11 disks (n = 18) for which |K(H_4)| > 2000.
- all_parts_summary.txt: raw summary line of each of the parts.
The lists of open disks (gzipped, 258 MB) are kept separately.

## Second stage (same day, 11:20-13:20): the open disks
- rerun_open.sh (in code/exhaustive_disks_C): the 23,935,868 disks not K-null within 5 steps were run again with 12 iterations
  and limit 3000 (rerun_summary.tsv). All became K-null except 28 (n = 18) that exceeded the limit.
- hard39_result.tsv: those 28 plus the 11 that exceeded 2000 in the first stage, run with limit 12000;
  hard2_result.tsv: the last 2, with limit 60000 (largest clique graph computed: 57,768 vertices).
- final_summary.tsv / final_table.tex: for every n, number of disks whose first 0-dismantlable iterate is K^i(B).

Conclusion: every flag disk B with at most 18 vertices (550,551,244 disks) is K-null; K^i(B) is 0-dismantlable
for some i <= 8, and every iterate K^n(B) (n >= 0) is 1-dismantlable (the cores of the iterates before the first
0-dismantlable one were certified 1-dismantlable; from that one on, all iterates are 0-dismantlable).
