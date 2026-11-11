#!/usr/bin/env bash
# Compila geng con el plug-in k1proc.c. Uso: bash compilar.sh <directorio_fuentes_nauty>
NS=${1:?directorio de nauty}
gcc -O3 -march=native -I"$NS" -DMAXN=WORDSIZE -DWORDSIZE=32 -DOUTPROC=k1out -DSUMMARY=k1sum \
  -o geng_k1 "$NS/geng.c" k1proc.c "$NS/gtools.c" "$NS/nauty.c" "$NS/nautil.c" "$NS/naugraph.c" "$NS/schreier.c" "$NS/naurng.c"
