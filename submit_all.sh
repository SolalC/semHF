#!/bin/bash
# Full pipeline per subset (female-vs-male subtraction, both orderings): 01 -> 02 -> 03 -> 04.
# Each subset writes to its own folder: out_base/{ALL,I,NI}.
# Usage: bash submit_all.sh            # ALL, I and NI
#        bash submit_all.sh I NI       # only these subsets
set -euo pipefail
cd /scratch/user/uqschauq/sem/Code
mkdir -p /scratch/user/uqschauq/sem/Code/logs
SUBSETS=("$@"); [ ${#SUBSETS[@]} -eq 0 ] && SUBSETS=(ALL I NI)

for s in "${SUBSETS[@]}"; do
  export HF_GSEM_SUBSET=$s
  j1=$(sbatch --parsable --export=ALL --job-name=hf_${s}_01_ldsc run_01.slurm)
  j2=$(sbatch --parsable --export=ALL --job-name=hf_${s}_02_userGWAS --dependency=afterok:${j1} run_02.slurm)
  j3=$(sbatch --parsable --export=ALL --job-name=hf_${s}_03_combine  --dependency=afterok:${j2} run_03.slurm)
  j4=$(sbatch --parsable --export=ALL --job-name=hf_${s}_04_miami    --dependency=afterok:${j3} run_04.slurm)
  echo "${s}: 01=${j1}  02=${j2} (array 1-22)  03=${j3}  04=${j4}"
done

