#!/bin/bash
#SBATCH --job-name=hf_gsem_01_ldsc
#SBATCH --account=a_imb_cpdg
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/user/uqschauq/sem/Code/logs/%x_%j.out
set -euo pipefail

SCRIPT_DIR=/scratch/user/uqschauq/sem/Code
module load r/4.4.2-heavy

export HF_GSEM_CONFIG="${SCRIPT_DIR}/00_config.R"
cd "${SCRIPT_DIR}"
Rscript 01_munge_ldsc_model.R
