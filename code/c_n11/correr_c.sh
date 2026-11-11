#!/usr/bin/env bash
# Búsqueda exhaustiva en 11 vértices con geng_k1 (geng + plug-in en C), en NUC procesos y 2000 partes.
# Reanudable; cada parte terminada se comprime y se respalda en $BK. Uso: bash correr_c.sh [NUC]
set -u; cd "$(dirname "$0")"
NUC=${1:-2}; N=11; PARTES=2000
DIR=$PWD/resultados_n11; BK=/mnt/user-data/outputs/n11c_respaldo; mkdir -p "$DIR" "$BK"
for f in "$BK"/parte_*.ok; do [ -e "$f" ] || continue; b=$(basename "$f" .ok)
  [ -e "$DIR/$b.ok" ] || { cp "$BK/$b.tsv.gz" "$DIR/" 2>/dev/null && cp "$f" "$DIR/"; }; done
rm -f "$DIR"/parte_*.log "$DIR"/parte_*.tsv
echo "inicio $(date '+%F %T')  partes listas: $(ls "$DIR"/*.ok 2>/dev/null | wc -l)/$PARTES"
export N PARTES DIR BK
seq 0 $((PARTES-1)) | xargs -P "$NUC" -I{} bash -c '
  i={}; f="$DIR/parte_$i"; [ -f "$f.ok" ] && exit 0
  if ./geng_k1 -c -q "$N" "$i/$PARTES" > "$f.tsv" 2> "$f.log"; then
    gzip -f "$f.tsv" && mv "$f.log" "$f.ok" && cp "$f.tsv.gz" "$f.ok" "$BK/"
    k=$(ls "$DIR"/*.ok | wc -l); [ $((k % 50)) -eq 0 ] && echo "$k/$PARTES partes  $(date "+%F %T")"
  fi'
echo "fin $(date '+%F %T')"
