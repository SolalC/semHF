#!/bin/bash
#SBATCH --job-name=hf_gsem_02_userGWAS
#SBATCH --account=a_imb_cpdg
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00                  # chr1 is the longest task; adjust after a first run
#SBATCH --array=1-22
#SBATCH --output=/scratch/user/uqschauq/sem/Code/logs/%x_%A_chr%a.out
set -euo pipefail

SCRIPT_DIR=/scratch/user/uqschauq/sem/Code
module load r/4.4.2-heavy

export HF_GSEM_CONFIG="${SCRIPT_DIR}/00_config.R"
cd "${SCRIPT_DIR}"
Rscript 02_userGWAS_chr.R
