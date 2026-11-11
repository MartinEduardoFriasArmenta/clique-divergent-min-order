"""Collects the parts written by run_exhaustive.sh and prints, for each order n, the number of flag disks, how many are
K-null with index i (first 0-dismantlable iterate K^i(B)), how many remain open after IT steps, exceed the limit,
or have some iterate whose core is not 1-dismantlable (NO) or undecided (?).  Writes exhaustive_summary.tsv and
exhaustive_table.tex.   Usage: python3 summarize.py [outdir]"""
import sys, glob, re, collections
d = sys.argv[1] if len(sys.argv) > 1 else 'out'
agg = collections.defaultdict(lambda: None); gen = collections.Counter(); parts = collections.Counter()
for f in glob.glob(d + '/n*_r*.sum'):
    n = int(re.search(r'n(\d+)_r', f).group(1))
    line = [l for l in open(f) if l.startswith('>S')][0].split()[1:]
    v = list(map(int, line))
    if agg[n] is None: agg[n] = v
    else:
        a = agg[n]; k = len(v)
        agg[n] = [a[i] + v[i] for i in range(k - 2)] + [max(a[k - 2], v[k - 2]), max(a[k - 1], v[k - 1])]
    parts[n] += 1
    g = open(f.replace('.sum', '.gen')).read()
    m = re.search(r'>Z (\d+)', g); gen[n] += int(m.group(1)) if m else 0
rows = []
for n in sorted(k for k in agg if agg[k]):
    v = agg[n]; IT = len(v) - 8
    total, knull, rest = v[0], v[1:IT + 2], v[IT + 2:]
    nopen, nbig, nno, nunk, maxK, maxC = rest
    assert total == gen[n], (n, total, gen[n])
    rows.append((n, parts[n], total, knull, nopen, nbig, nno, nunk, maxK, maxC))
with open('exhaustive_summary.tsv', 'w') as f:
    f.write('n\tparts\tdisks\t' + '\t'.join('Knull_%d' % i for i in range(len(rows[0][3]))) + '\topen\t>L\tNO\t?\tmax|K(H)|\tmax|core|\n')
    for r in rows: f.write('\t'.join(map(str, [r[0], r[1], r[2], *r[3], *r[4:]])) + '\n')
with open('exhaustive_table.tex', 'w') as f:
    k = len(rows[0][3])
    f.write('\\begin{tabular}{r r ' + 'r' * k + ' rrr}\\toprule\n$n$ & disks & ' + ' & '.join('$K$-null, index %d' % i for i in range(k))
            + ' & open & $>L$ & not 1-dism.\\\\\\midrule\n')
    for r in rows: f.write(' & '.join(map(str, [r[0], r[2], *r[3], r[4], r[5], r[6] + r[7]])) + '\\\\\n')
    f.write('\\bottomrule\\end{tabular}\n')
for r in rows: print(r)
