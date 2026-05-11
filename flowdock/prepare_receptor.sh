#!/bin/bash

INPUT_DIR="PDB"
OUTPUT_DIR="clean_PDB"

mkdir -p "$OUTPUT_DIR"

while read pdb lig chain; do

infile="$INPUT_DIR/${pdb}.pdb"

receptor_out="$OUTPUT_DIR/${pdb}_receptor.pdb"
ligand_pdb="$OUTPUT_DIR/${pdb}_ligand.pdb"
ligand_sdf="$OUTPUT_DIR/${pdb}_ligand.sdf"

echo "Procesando $pdb | Ligando=$lig | Cadena=$chain"

awk -v CHAIN="$chain" '

/^ATOM/ {
    if (substr($0,22,1)==CHAIN) print
    next
}

/^TER/ {print}
/^END/ {print}

' "$infile" > "$receptor_out"

awk -v LIG="$lig" -v CHAIN="$chain" '

/^HETATM/ {
    resname=substr($0,18,3)
    ch=substr($0,22,1)

    if (resname==LIG && ch==CHAIN) print
    next
}

/^CONECT/ {print}
/^END/ {print}

' "$infile" > "$ligand_pdb"

obabel "$ligand_pdb" -O "$ligand_sdf" >/dev/null 2>&1

echo "  -> receptor creado"
echo "  -> ligando creado"
echo ""

done < ligands.txt

