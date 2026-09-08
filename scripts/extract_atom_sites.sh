#!/usr/bin/env bash

# usage
# extract_atom_sites.sh bulk_afdb_path.tar organism_name

BULK_AFDB_FILE="$1" # e.g., data/afdb_tar/XXX.tar
ORGANISM="$2" # e.g., fruit_fly

# Create directories
EXTRACTED_DIR="data/afdb_extracted/$ORGANISM/"
ATOM_SITE_DIR="data/atom_sites/$ORGANISM/"
mkdir -p "$EXTRACTED_DIR" "$ATOM_SITE_DIR"

# Extract .tar archive into EXTRACTED_DIR
tar -xf "$BULK_AFDB_FILE" -C "$EXTRACTED_DIR"

# Unzip all .pdb.gz files
gzip -d "$EXTRACTED_DIR"*.pdb.gz

# Extract atomsite data from each .pdb file
for FILE in "$EXTRACTED_DIR"*.pdb; do

    OUT_CSV="$ATOM_SITE_DIR/$(basename "$FILE" .pdb).csv"
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
