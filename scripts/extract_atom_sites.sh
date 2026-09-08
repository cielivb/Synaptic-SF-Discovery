#!/usr/bin/env bash

# usage
# extract_atom_sites.sh bulk_afdb_path.tar atom_site_outdir

BULK_AFDB_FILE="$1"
STRUCTURE_DIR="$2"

# Create directories
PDB_DIR="$STRUCTURE_DIR/pdb"
ATOM_SITE_DIR="$STRUCTURE_DIR/atom_site"
mkdir -p "$PDB_DIR" "$ATOM_SITE_DIR"

# Extract .tar archive into PDB_DIR
tar -xf "$BULK_AFDB_FILE" -C "$PDB_DIR"

# Remove .cif.gz files
rm -f "$PDB_DIR"/*.cif.gz

# Unzip all .pdb.gz files
gzip -d "$PDB_DIR"/*.pdb.gz

# Extract atomsite data from each .pdb file
for FILE in "$PDB_DIR"/*.pdb; do

    OUT_CSV="$ATOM_SITE_DIR/$(basename "$FILE" .pdb).csv"
    echo "residue,x,y,z,pLDDT" > "$OUT_CSV" # header line

    # PDB coordinates lines start with ATOM and the columns are whitespace-
    # separated. Fields 7, 8, 9 correspond to x, y, z coords of an atom.
    # Field 11 corresponds to the pLDDT score of the residue in AFDB PDB files.
    # Field 6 is the residue ID.
    awk '$1=="ATOM"{printf("%s,%s,%s,%s,%s\n", $6, $7, $8, $9, $11)}' "$FILE" >> "$OUT_CSV"

done
