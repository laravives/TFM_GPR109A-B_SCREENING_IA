#!/bin/bash -l

#SBATCH --job-name=Amber
#SBATCH --account=p201228
#SBATCH --partition=gpu
#SBATCH --qos=default
#SBATCH --time=24:00:00

module use /project/home/p201228/install_lxp/2025.1/modules/all/
module load Amber/24.3-foss-2025a-CUDA-12.8.0-pmemd

source /project/home/p201228/miniconda3/etc/profile.d/conda.sh
conda activate amber

while read -r base; do

cd $base

pmemd.cuda -O -i min1.in -p complex.prmtop -c complex.inpcrd -ref complex.inpcrd -o min1.out -r min1.rst
pmemd.cuda -O -i min2.in -p complex.prmtop -c complex.inpcrd -ref complex.inpcrd -o min2.out -r min2.rst
pmemd.cuda -O -i heat.in -p complex.prmtop -c min2.rst -ref min2.rst -o heat.out -r heat.rst -x heat.nc
pmemd.cuda -O -i eq.in -p complex.prmtop -c heat.rst -ref heat.rst -o eq.out -r eq.rst -x eq.nc
pmemd.cuda -O -i prod.in -p complex.prmtop -c eq.rst -o prod.out -r prod.rst -x prod.nc

cpptraj complex.prmtop <<EOF
trajin prod.nc
strip :WAT,Na+,Cl-
trajout prod_dry.nc netcdf
EOF

cpptraj complex.prmtop <<EOF
parm complex.prmtop
parmstrip :WAT,Na+,Cl-
parmbox nobox
parmwrite out complex_dry.prmtop
EOF


cat > rmsd.in <<EOF
parm complex.prmtop
trajin prod.nc
reference complex.inpcrd

rms BackboneRMSD reference @CA,C,N out rmsd_backbone.dat

run
quit
EOF

cpptraj -i rmsd.in

cd ..

done < lig.txt
