#!/bin/bash
#SBATCH --job-name=hf_gsem_03_combine
#SBATCH --account=a_XXXX                 # EDIT: your Bunya account
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/user/uqschauq/sem/hf_sex_gwas_by_subtraction/logs/%x_%j.out
set -euo pipefail

SCRIPT_DIR=/scratch/user/uqschauq/sem/Code
module load r/4.4.0                              # EDIT

export HF_GSEM_CONFIG="${SCRIPT_DIR}/00_config.R"
cd "${SCRIPT_DIR}"
Rscript 03_combine_qc.R
