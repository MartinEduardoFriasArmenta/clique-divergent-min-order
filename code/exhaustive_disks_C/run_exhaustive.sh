#!/bin/bash
# Exhaustive experiment of Section 7: all flag disks with NMIN <= n <= NMAX vertices.
# Usage: ./run_exhaustive.sh NMIN NMAX [JOBS=number of cores] [MOD=parts per order n] [IT=5] [LIMIT=2000]
# Each part is  gen_disks_ca n r MOD | kdisk -S -O -i IT -L LIMIT  and is resumable: a part whose out/n<n>_r<r>.sum
# already exists is skipped, so the script can be stopped and relaunched at any time.
# The exceptional and "open" disks of each part go (gzipped) to out/n<n>_r<r>.tsv.gz.
# Then:  python3 summarize.py out   (tables in TSV and LaTeX).
set -u
NMIN=$1; NMAX=$2; JOBS=${3:-$(nproc)}; MOD=${4:-64}; IT=${5:-5}; LIMIT=${6:-2000}
DIR=$(cd "$(dirname "$0")" && pwd); mkdir -p out
run_part() {   # n r
  local n=$1 r=$2 mod=$MOD
  [ "$n" -le 10 ] && mod=1
  [ "$r" -ge "$mod" ] && return 0
  local base=out/n${n}_r${r}
  [ -s "$base.sum" ] && return 0
  "$DIR/gen_disks_ca" "$n" "$r" "$mod" 2> "$base.gen" | "$DIR/kdisk" -S -O -i "$IT" -L "$LIMIT" 2> "$base.tmp" | gzip > "$base.tsv.gz"
  grep -q '^>S' "$base.tmp" && mv "$base.tmp" "$base.sum"
}
export -f run_part; export DIR MOD IT LIMIT
for n in $(seq "$NMIN" "$NMAX"); do for r in $(seq 0 $((MOD - 1))); do echo "$n $r"; done; done \
  | xargs -P "$JOBS" -n 2 bash -c 'run_part "$0" "$1"'
echo "done: $(ls out/*.sum 2>/dev/null | wc -l) parts finished"
