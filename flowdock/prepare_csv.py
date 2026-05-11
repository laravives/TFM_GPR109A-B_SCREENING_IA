#!/usr/bin/env python3

import csv
import os
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

# CONFIGURACION
BASE_DIR = Path("/project/home/p201228/FlowDock/GPR109A")
INPUT_CSV = BASE_DIR / "inputsA.csv"
OUTPUT_CSV = BASE_DIR / "inputs.csv"
PDB_DIR = BASE_DIR / "clean_PDB"
SDF_DIR = BASE_DIR / "ligands_sdf"
SDF_DIR.mkdir(exist_ok=True)

# AA MAPA
aa3to1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
    "MSE":"M"
}

# EXTRAER SECUENCIA PDB
def pdb_to_sequence(pdb_file):
    seq = []
    seen = set()

    with open(pdb_file) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue

            resname = line[17:20].strip()
            chain = line[21].strip()
            resnum = line[22:26].strip()

            key = (chain, resnum)
            if key in seen:
                continue
            seen.add(key)

            aa = aa3to1.get(resname, "X")
            seq.append(aa)

    return "".join(seq)

# SMILES -> SDF
def smiles_to_sdf(smiles, outfile):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    mol = Chem.AddHs(mol)

    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)

    writer = Chem.SDWriter(str(outfile))
    writer.write(mol)
    writer.close()

# PROCESAR
rows_out = []

with open(INPUT_CSV) as f:
    reader = csv.DictReader(f)

    for row in reader:
        ID = row["id"].strip()
        smiles = row["input_ligand"].strip()

        pdb_file = PDB_DIR / f"{ID}_receptor.pdb"
        if not pdb_file.exists():
            print(f"[WARNING] Missing PDB for {ID}: {pdb_file}")
            continue

        # sequence
        seq = pdb_to_sequence(pdb_file)

        # sdf ligand
        sdf_file = SDF_DIR / f"lig_{ID}.sdf"
        smiles_to_sdf(smiles, sdf_file)

        rows_out.append({
            "id": ID,
            "input_receptor": seq,
            "input_ligand": str(sdf_file.relative_to(BASE_DIR)),
            "input_template": str(pdb_file.relative_to(BASE_DIR))
        })

        print(f"[OK] {ID}")

# ESCRIBIR OUTPUT CSV
with open(OUTPUT_CSV, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "id",
            "input_receptor",
            "input_ligand",
            "input_template"
        ]
    )
    writer.writeheader()
    writer.writerows(rows_out)


