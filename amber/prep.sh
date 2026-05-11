#!/bin/bash

source /project/home/p201228/miniconda3/etc/profile.d/conda.sh
conda activate amber

# Leer todas las carpetas desde lig.txt
while read dir; do

    echo "Procesando $dir ..."

    cd "$dir" || continue

    # Limpiar con pdb4amber

    pdb4amber -i $dir.pdb -o complex.pdb --reduce

    # Separar proteína y ligando

    grep "^ATOM" complex.pdb > protein.pdb
    grep "^HETATM" complex.pdb > ligand.pdb

    # Preparar ligando

    antechamber \
        -i ligand.pdb \
        -fi pdb \
        -o ligand.mol2 \
        -fo mol2 \
        -c bcc \
        -at gaff2 \
        -s 2 \
        -nc 0

    # Generar frcmod

    parmchk2 \
        -i ligand.mol2 \
        -f mol2 \
        -o ligand.frcmod

    # Crear tleap.in

    cat > tleap.in <<EOF
source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p

LIG = loadmol2 ligand.mol2
loadamberparams ligand.frcmod

PROT = loadpdb protein.pdb

COMPLEX = combine {PROT LIG}

check COMPLEX

solvatebox COMPLEX TIP3PBOX 10.0
addions COMPLEX Cl- 0

saveamberparm COMPLEX complex.prmtop complex.inpcrd
saveamberparm PROT receptor.prmtop receptor.inpcrd
saveamberparm LIG ligand.prmtop ligand.inpcrd

quit
EOF

    # Ejecutar tleap

    tleap -f tleap.in > tleap.log

    # 8. Volver atrás

    cd ..

done < lig.txt
