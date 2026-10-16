# semHF: sex-specific heart failure GWAS by subtraction (Genomic SEM)

Cholesky model on the female (HF_F) and male (HF_M) all-cause heart failure GWAS:

    g1 = l11 F1
    g2 = l21 F1 + l22 F2        Var(F1) = Var(F2) = 1, Cov(F1, F2) = 0

F2 is the genetic component of trait 2 not shared with trait 1. Both orderings are run:
`F_specific` (male subtracted, F2 = female-specific) and `M_specific` (female subtracted,
F2 = male-specific).

## Pipeline

| Step | Script | What it does |
|---|---|---|
| 00 | `00_config.R` | Paths and settings shared by all steps (edit the marked lines) |
| 01 | `01_munge_ldsc_model.R` | Effective N, munge, LDSC, Cholesky fit and power summary, `sumstats()` |
| 01b | `01b_test_specific_factor.R` | Optional: chi-square test of r_g = 1 (one factor vs Cholesky) |
| 02 | `02_userGWAS_chr.R` | `userGWAS` per chromosome (SLURM array 1-22), F1 ~ SNP and F2 ~ SNP |
| 03 | `03_combine_qc.R` | Combine chromosomes, repair multi-line records, drop errors and smoothed SNPs, QC, extract SNPs below `sig_alpha` |
| 04 | `04_miami_plot.R` | Miami plot (F1 top, F2 bottom), peaks marked (filled: p < 5e-8; hollow: p < `sig_alpha`) |

Run everything with `bash submit_all.sh` (01 -> 02 -> 03 -> 04, SLURM dependencies).
Before the first run, set the reference paths in `00_config.R` and the account and R module in `run_0*.slurm`.

## Main outputs (in `out_dir`)

- `ldsc_summary.tsv`, `cholesky_summary.tsv`: h2, r_g, loadings, specific h2, power
- `{ordering}_{F1,F2}_genomewide.tsv.gz`: all SNPs; `pass` = no error and S not smoothed
- `{ordering}_{F1,F2}_significant_p1e-06.tsv`: passing SNPs with p < `sig_alpha`
- `userGWAS_qc.tsv`: per ordering x factor counts, lambda_GC, smoothed SNPs dropped, z check against the closed form
- `plots/{ordering}_miami.png|pdf`, `plots/{ordering}_miami_peaks.tsv`
