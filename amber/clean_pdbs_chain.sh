#!/bin/bash

INPUT_DIR="PDB"
OUTPUT_DIR="clean_PDB"

mkdir -p $OUTPUT_DIR

while read pdb lig chain; do
    infile="$INPUT_DIR/${pdb}.pdb"
    outfile="$OUTPUT_DIR/${pdb}_clean.pdb"

    echo "Procesando $pdb (ligando: $lig, cadena: $chain)"

    awk -v LIG="$lig" -v CHAIN="$chain" '

    /^ATOM/ {
        if (substr($0,22,1) == CHAIN) {
            print
        }
        next
    }

    /^HETATM/ {
        resname=substr($0,18,3)
        ch=substr($0,22,1)

        if (resname == LIG && ch == CHAIN) {
            print
        }
        next
    }

    /^TER/ {print}
    /^END/ {print}

    ' "$infile" > "$outfile"

done < ligands.txt

echo "Limpieza completa terminada"
