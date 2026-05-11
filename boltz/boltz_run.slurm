#!/bin/bash
#SBATCH --job-name=boltz_run
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=04:00:00
#SBATCH --account=p201228
#SBATCH --qos=short
#SBATCH --output=boltz_run_%j.out

source /project/home/p201228/miniconda3/etc/profile.d/conda.sh
conda activate boltz_env

cd /project/home/p201228/boltz_runs/GPR109A_B_ligandos

for yaml in *.yaml; do
    boltz predict $yaml \
        --model boltz2 \
        --devices 1 \
        --use_msa_server \
        --affinity_mw_correction \
        --override
done
