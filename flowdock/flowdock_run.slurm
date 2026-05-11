#!/bin/bash
#SBATCH -A p201228
#SBATCH -p gpu
#SBATCH -q default
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --array=1-20
#SBATCH -J flowdock
#SBATCH -o logs/flowdock_%A_%a.out
#SBATCH -e logs/flowdock_%A_%a.err

cd /project/home/p201228/FlowDock/GPR109A
mkdir -p logs results

source /etc/profile
module purge
module load env/release/2025.1
module load Apptainer/1.4.2-GCCcore-14.2.0

LINE=$((SLURM_ARRAY_TASK_ID))

readarray -t FIELDS < <(
python3 - <<PY
import csv
with open("inputs.csv") as f:
    rows=list(csv.reader(f))
row=rows[$LINE]
for x in row:
    print(x)
PY
)

ID="${FIELDS[0]}"
RECEPTOR="${FIELDS[1]}"
LIGAND="${FIELDS[2]}"
TEMPLATE="${FIELDS[3]}"

OUTDIR=/project/home/p201228/FlowDock/GPR109A/results/$ID
mkdir -p "$OUTDIR"

echo "Running task $SLURM_ARRAY_TASK_ID"
echo "ID=$ID"
echo "Ligand=$LIGAND"

apptainer exec --nv \
-B /project/home/p201228:/project/home/p201228 \
/project/home/p201228/flowdock_latest.sif bash -lc "

source /opt/conda/etc/profile.d/conda.sh
conda activate flowdock

python /project/home/p201228/FlowDock/flowdock/sample.py \
ckpt_path=/project/home/p201228/FlowDock/checkpoints/checkpoints/esmfold_prior_paper_weights-EMA.ckpt \
model.cfg.prior_type=esmfold \
sampling_task=batched_structure_sampling \
sample_id='$ID' \
input_receptor='$RECEPTOR' \
input_ligand='$LIGAND' \
input_template='$TEMPLATE' \
out_path='$OUTDIR' \
n_samples=5 \
num_steps=20 \
chunk_size=1 \
trainer=gpu
"
