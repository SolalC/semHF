#!/bin/bash
# Full pipeline, both orderings (F_specific, M_specific): 01 -> 02 -> 03 -> 04.
set -euo pipefail
cd /scratch/user/uqschauq/sem/Code
mkdir -p /scratch/user/uqschauq/sem/hf_sex_gwas_by_subtraction/logs
j1=$(sbatch --parsable run_01.slurm)
j2=$(sbatch --parsable --dependency=afterok:${j1} run_02.slurm)
j3=$(sbatch --parsable --dependency=afterok:${j2} run_03.slurm)
j4=$(sbatch --parsable --dependency=afterok:${j3} run_04.slurm)
echo "Submitted: 01=${j1}  02=${j2} (array 1-22)  03=${j3}  04=${j4}"
