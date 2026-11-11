#!/bin/bash
# Re-runs the "open" disks (not K-null within 5 steps, limit not exceeded) of each part with IT iterations.
# Usage: ./rerun_open.sh [JOBS=2] [IT=12] [LIMIT=3000]
cd "$(dirname "$0")"; JOBS=${1:-2}; IT=${2:-12}; LIMIT=${3:-3000}
part() { f=$1; b=out2/$(basename "$f" .tsv.gz); [ -s "$b.sum" ] && return 0
  zcat "$f" | awk -F'\t' '{k=0; for(i=4;i<=NF;i++){if($i=="K1"||substr($i,1,1)==">")k=1}; if(!k) print $1}' > "$b.in"
  ./kdisk -S -O -i "$IT" -L "$LIMIT" < "$b.in" 2> "$b.tmp" | gzip > "$b.tsv.gz"
  grep -q '^>S' "$b.tmp" && mv "$b.tmp" "$b.sum" && rm -f "$b.in"; }
export -f part; export IT LIMIT
ls out/n*_r*.tsv.gz | sort -V | xargs -P "$JOBS" -n 1 bash -c 'part "$0"'
echo FIN > fin2.txt
