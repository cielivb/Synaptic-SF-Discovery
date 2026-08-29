#!/usr/bin/env bash

# usage
# extract_coordinates.sh data_pdb data_csv

DATA_PDB_DIR=$1
DATA_CSV_DIR=$2


echo "$(date +%H:%M:%S) Extracting coordinates from .pdb files into .csv files ..."

PDB_FILES="$DATA_PDB_DIR/UP000005640/"*.pdb
for FILE in PDB_FILES; do

    OUT_CSV="$DATA_CSV_DIR/"$(basename "$FILE" .pdb).csv
    echo "x,y,z" > "$OUT_CSV"

    # PDB coordinates lines start with ATOM and the columns are whitespace-
    # separated. Fields 7, 8, 9 correspond to x, y, z coords of an atom.
    awk '$1=="ATOM"{printf("%s,%s,%s\n", $7, $8, $9)}' "$PDB" >> "$OUT_CSV"

done

echo "$(date +%H:%M:%S) Coordinates extracted"
