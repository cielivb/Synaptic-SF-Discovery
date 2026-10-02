#!/usr/bin/env bash

# usage
# extract_atom_sites.sh pdb_dir csv_dir

PDB_DIR=$1
CSV_DIR=$2

# Extract atomsite data from each .pdb file in PDB_DIR
for FILE in "$PDB_DIR"/*.pdb; do

    UNIPROT=$(basename "$FILE" | cut -d'-' -f2)

    OUT_CSV="$CSV_DIR/$UNIPROT.csv"
    echo "residue,x,y,z,pLDDT" > "$OUT_CSV" # header line

    awk '
    substr($0,1,6)=="ATOM  " {
        residue = substr($0,23,4) + 0
        x   = substr($0,31,8) + 0
        y   = substr($0,39,8) + 0
        z   = substr($0,47,8) + 0
        plddt = substr($0,61,6) + 0
        printf("%d,%f,%f,%f,%f\n", residue, x, y, z, plddt)
    }
    ' "$FILE" >> "$OUT_CSV"

done
