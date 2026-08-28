#!/usr/bin/env bash

# Usage:
# atomsite2csv.sh input.pdb output.csv

PDB="$1"
OUT_CSV="$2"

echo "x,y,z" > "$OUT_CSV"

# PDB coordinates lines start with ATOM and the columns are whitespace-
# separated. Fields 7, 8, 9 correspond to x, y, z coords of an atom.
awk '$1=="ATOM"{printf("%s,%s,%s\n", $7, $8, $9)}' "$PDB" >> "$OUT_CSV"
