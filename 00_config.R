## =============================================================================
## GWAS-by-subtraction (Genomic SEM Cholesky model): sex-stratified heart failure
## 00_config.R -- settings shared by all scripts. EDIT THE MARKED VARIABLES.
## =============================================================================

## ---- USER SETTINGS ---------------------------------------------------------
sumstats_dir <- "/QRISdata/Q4486/gwas_cvd_sumstats/HF_sex_stratified_CJ/20260609"
out_base     <- "/scratch/user/uqschauq/sem/hf_sex_gwas_by_subtraction"   # one subfolder per subset

## HF subsets. Each runs the same female-vs-male subtraction on
## gcta_{subset}_HF_EU_FEMALE and gcta_{subset}_HF_EU_MALE. The subset for a run comes
## from the environment (HF_GSEM_SUBSET, set by submit_all.sh); default ALL.
subsets <- c("ALL", "I", "NI")
subset  <- Sys.getenv("HF_GSEM_SUBSET", "ALL"); if (!nzchar(subset)) subset <- "ALL"
stopifnot(subset %in% subsets)
out_dir     <- file.path(out_base, subset)
file_female <- file.path(sumstats_dir, paste0("gcta_", subset, "_HF_EU_FEMALE"))
file_male   <- file.path(sumstats_dir, paste0("gcta_", subset, "_HF_EU_MALE"))

## Genomic SEM / LDSC reference files (EDIT)
hm3_snplist  <- "/QRISdata/Q9564/LD/LDSC/GRCh37/eur_w_ld_chr/w_hm3.snplist"
ld_dir       <- "/QRISdata/Q9564/LD/LDSC/GRCh37/eur_w_ld_chr/"   # keep trailing slash
ref_1000g    <- "/QRISdata/Q9564/LD/1000G/GRCh37/reference.1000G.maf.0.005.txt.gz"   # gz read by data.table::fread (needs R.utils)

## Case/control counts per subset and sex. Leave NA to estimate the effective N from
## the SEs (Neff = 4 / (2 p q se^2), median over common SNPs; ignores covariates).
counts <- list(            # c(ncase_female, ncontrol_female, ncase_male, ncontrol_male)
  ALL = c(NA, NA, NA, NA),
  I   = c(NA, NA, NA, NA),
  NI  = c(NA, NA, NA, NA)
)
ncase_female <- counts[[subset]][1]; ncontrol_female <- counts[[subset]][2]
ncase_male   <- counts[[subset]][3]; ncontrol_male   <- counts[[subset]][4]

## Population prevalence per subset and sex, for the liability-scale conversion.
## Leave NA to work on the observed scale. h2, loadings, b and the F2 odds ratios
## depend on this choice; with different prevalences by sex, so do the F2 z-statistics.
pop_prev <- list(          # c(female, male)
  ALL = c(NA, NA),
  I   = c(NA, NA),
  NI  = c(NA, NA)
)
pop_prev_female <- pop_prev[[subset]][1]
pop_prev_male   <- pop_prev[[subset]][2]

maf_filter  <- 0.01
sig_alpha   <- 1e-6     # threshold for extracting SNPs (1e-6 = suggestive; 5e-8 = genome-wide)
z_smooth_max <- 0       # drop SNPs whose S matrix had to be smoothed (Z_smooth > this); 0 = drop any smoothed SNP
info_filter <- 0.9      # no INFO column in these files; munge skips it
seed        <- 20261015
## ---------------------------------------------------------------------------

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
traits  <- c("HF_F", "HF_M")

## Cholesky orderings. t1 = trait subtracted; t2 = trait whose residual factor F2
## is the sex-specific component.
orderings <- list(
  F_specific = c(t1 = "HF_M", t2 = "HF_F"),   # female signal not shared with males
  M_specific = c(t1 = "HF_F", t2 = "HF_M")    # male signal not shared with females
)


## Cholesky model: g1 = l11 F1 ; g2 = l21 F1 + l22 F2
cholesky_model <- function(t1, t2, with_snp = FALSE) {
  m <- stringr::str_glue(
    "F1 =~ NA*{t1} + {t2}
     F2 =~ NA*{t2}
     F1 ~~ 1*F1
     F2 ~~ 1*F2
     F1 ~~ 0*F2
     {t1} ~~ 0*{t1}
     {t2} ~~ 0*{t2}
     {t1} ~~ 0*{t2}")
  if (with_snp) m <- paste(m, "F1 ~ SNP", "F2 ~ SNP", sep = "\n")
  as.character(m)
}
