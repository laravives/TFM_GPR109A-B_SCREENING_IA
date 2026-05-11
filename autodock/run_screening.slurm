#!/bin/bash
#SBATCH --job-name=adgpu
#SBATCH --account=p201228
#SBATCH --partition=gpu
#SBATCH --qos=default
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gpus-per-node=4
#SBATCH --time=02:00:00

source "/project/home/p201228/miniconda3/etc/profile.d/conda.sh"
conda activate autodock

BASE_DIR="/project/home/p201228/gpcr_screening"
ADGPU="$BASE_DIR/bin/autodock_gpu"
INDEX_FILE="$BASE_DIR/analysis/receptor_index.txt"

ENTRY=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$INDEX_FILE")
PROTEIN="${ENTRY%% *}"
PDB_ID="${ENTRY##* }"

declare -A CHAIN_MAP
CHAIN_MAP["8J6P"]="R"
CHAIN_MAP["8J6Q"]="R"
CHAIN_MAP["8JZ7"]="A"
CHAIN_MAP["9IZC"]="A"
CHAIN_MAP["8IJA"]="A"
CHAIN_MAP["8J6J"]="R"
CHAIN_MAP["8JEF"]="A"
CHAIN_MAP["8JEI"]="A"
CHAIN_MAP["9JID"]="A"
CHAIN=${CHAIN_MAP[$PDB_ID]}

MAPS_DIR="$BASE_DIR/$PROTEIN/receptors/maps_blind_chain"
FLD="$MAPS_DIR/${PDB_ID}_blind_chain${CHAIN}.maps.fld"
LIG_DIR="$BASE_DIR/common"
OUT_DIR="$BASE_DIR/$PROTEIN/results_blind_chain"
mkdir -p "$OUT_DIR/poses" "$OUT_DIR/dlg" "$OUT_DIR/logs"

[[ ! -f "$FLD" ]] && { echo "ERROR: no hay $FLD"; exit 1; }

nvidia-smi --query-gpu=index,name --format=csv,noheader
N_GPUS=$(nvidia-smi --query-gpu=name --format=csv,noheader | wc -l)

LIGANDS=("$LIG_DIR"/*.pdbqt)
N_LIGS=${#LIGANDS[@]}

GPU_ID=0
PIDS=()
for lig_file in "${LIGANDS[@]}"; do
    lig_name=$(basename "$lig_file" .pdbqt)
    out_prefix="$OUT_DIR/poses/${PDB_ID}_${lig_name}"
    log_file="$OUT_DIR/logs/${PDB_ID}_${lig_name}.log"

    if [[ -f "${out_prefix}.xml" ]] && [[ -s "${out_prefix}.xml" ]]; then
        echo "03_run_screening_ciego.slurm"
        continue
    fi

    echo "  [GPU$GPU_ID] $PDB_ID + $lig_name"
    CUDA_VISIBLE_DEVICES=$GPU_ID "$ADGPU" \
        --ffile "$FLD" \
        --lfile "$lig_file" \
        --resnam "$out_prefix" \
        --nrun 50 \
        --lsmet ad \
        --heuristics 1 \
        --nev 2000000 \
        --gbest 1 \
        > "$log_file" 2>&1 &

    PIDS+=($!)
    GPU_ID=$(( (GPU_ID + 1) % N_GPUS ))

    # Máximo 4 jobs paralelos
    while [[ $(jobs -r | wc -l) -ge $N_GPUS ]]; do sleep 3; done
done

# Esperar todos
for pid in "${PIDS[@]}"; do wait "$pid"; done

# Mover .dlg a su directorio
for dlg in "$OUT_DIR/poses/${PDB_ID}_"*.dlg; do
    [[ -f "$dlg" ]] && mv "$dlg" "$OUT_DIR/dlg/"
done

