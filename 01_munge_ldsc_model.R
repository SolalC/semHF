## 01: harmonise sumstats, munge, LDSC, fit the Cholesky model (both orderings),
##     report paths / specific h2 / power, and prepare SNP files for userGWAS.
suppressPackageStartupMessages({ library(GenomicSEM); library(tidyverse) })
source(Sys.getenv("HF_GSEM_CONFIG", "00_config.R"))
set.seed(seed)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)

## ---- 1. Read GCTA-format files, effective N, standard column names ----------
read_gcta <- function(path) {
  read_table(path, col_types = cols(SNP = "c", A1 = "c", A2 = "c", .default = "d"))
}
neff_from_se <- function(d) {
  d %>%
    filter(freq > 0.1, freq < 0.9, is.finite(se), se > 0) %>%
    mutate(neff = 4 / (2 * freq * (1 - freq) * se^2)) %>%
    pull(neff) %>%
    median()
}
neff_from_counts <- function(ncase, ncontrol) 4 / (1 / ncase + 1 / ncontrol)

raw <- list(HF_F = read_gcta(file_female), HF_M = read_gcta(file_male))

neff <- c(
  HF_F = if (is.na(ncase_female)) neff_from_se(raw$HF_F) else neff_from_counts(ncase_female, ncontrol_female),
  HF_M = if (is.na(ncase_male))   neff_from_se(raw$HF_M) else neff_from_counts(ncase_male, ncontrol_male)
)
message("Effective N used: ", paste(names(neff), round(neff), sep = " = ", collapse = "; "))

clean_files <- c(HF_F = "HF_F_clean.txt", HF_M = "HF_M_clean.txt")
iwalk(raw, function(d, tr) {
  d %>%
    transmute(SNP, A1, A2, MAF = pmin(freq, 1 - freq), BETA = b, SE = se, P = p, N = neff[[tr]]) %>%
    filter(MAF >= maf_filter, is.finite(BETA), is.finite(SE), SE > 0) %>%
    write_tsv(clean_files[[tr]])
})
rm(raw); gc()

## ---- 2. Munge and LDSC ------------------------------------------------------
munge(files = unname(clean_files), hm3 = hm3_snplist, trait.names = traits,
      N = unname(neff), info.filter = info_filter, maf.filter = maf_filter)

liability   <- !is.na(pop_prev_female) && !is.na(pop_prev_male)
sample_prev <- if (liability) c(0.5, 0.5) else c(NA, NA)   # 0.5 because N = Neff
pop_prev    <- if (liability) c(pop_prev_female, pop_prev_male) else c(NA, NA)

covstruc <- ldsc(traits = paste0(traits, ".sumstats.gz"),
                 sample.prev = sample_prev, population.prev = pop_prev,
                 ld = ld_dir, wld = ld_dir, trait.names = traits,
                 stand = TRUE, ldsc.log = "HF_sex_ldsc")
saveRDS(covstruc, "covstruc.rds")

## LDSC summaries. V is ordered as vech(S): S11, S21, S22.
S <- covstruc$S; V <- covstruc$V
ldsc_summary <- tibble(
  scale    = if (liability) "liability" else "observed",
  h2_F     = S[1, 1], h2_F_se = sqrt(V[1, 1]),
  h2_M     = S[2, 2], h2_M_se = sqrt(V[3, 3]),
  rg       = covstruc$S_Stand[2, 1], rg_se = sqrt(covstruc$V_Stand[2, 2]),
  cross_trait_intercept = covstruc$I[2, 1]   # expected ~0: no sample overlap between sexes
) %>%
  mutate(z_rg_eq_1 = (1 - rg) / rg_se,
         p_rg_eq_1 = 2 * pnorm(-abs(z_rg_eq_1)),
         rg_l95 = rg - 1.96 * rg_se, rg_u95 = rg + 1.96 * rg_se)
print(ldsc_summary)

## ---- 3. Cholesky model, both orderings ---------------------------------------
idx <- set_names(seq_along(traits), traits)

summarise_ordering <- function(ord, label) {
  t1 <- ord[["t1"]]; t2 <- ord[["t2"]]
  fit <- usermodel(covstruc, estimation = "DWLS",
                   model = cholesky_model(t1, t2), CFIcalc = FALSE)
  loadings <- fit$results %>% filter(op == "=~")
  get_l <- function(f, t) loadings %>% filter(lhs == f, rhs == t)
  fit_l11 <- get_l("F1", t1); fit_l21 <- get_l("F1", t2); fit_l22 <- get_l("F2", t2)

  ## Closed-form check: l11 = sqrt(S11), l21 = S21 / l11, l22 = sqrt(S22 - l21^2)
  i <- idx[[t1]]; j <- idx[[t2]]
  s11 <- S[i, i]; s21 <- S[i, j]; s22 <- S[j, j]
  l11_cf <- sqrt(s11); l21_cf <- s21 / l11_cf; l22_cf <- sqrt(s22 - l21_cf^2)

  ## Specific h2 = S22 - S21^2 / S11, SE by delta method with the full V
  ## (includes the covariances between h2 and genetic covariance estimates).
  vech_pos <- function(a, b) if (a == b) c(1, 3)[a] else 2   # vech(S) = S11, S21, S22
  pos  <- c(vech_pos(i, i), vech_pos(i, j), vech_pos(j, j))   # S11, S21, S22 of this ordering
  grad <- c(s21^2 / s11^2, -2 * s21 / s11, 1)
  h2_spec    <- s22 - s21^2 / s11
  h2_spec_se <- sqrt(drop(t(grad) %*% V[pos, pos] %*% grad))

  ## Power of the subtraction GWAS (no sample overlap between sexes).
  b  <- l21_cf / l11_cf
  n1 <- neff[[t1]]; n2 <- neff[[t2]]
  rg <- ldsc_summary$rg
  ncp <- function(alpha, power = 0.8) (qnorm(1 - alpha / 2) + qnorm(power))^2
  sigma_u2 <- 1 / n2 + b^2 / n1

  tibble(
    ordering = label, subtracted = t1, residual_trait = t2,
    l11 = fit_l11$Unstand_Est, l11_se = fit_l11$Unstand_SE,
    l21 = fit_l21$Unstand_Est, l21_se = fit_l21$Unstand_SE,
    l22 = fit_l22$Unstand_Est, l22_se = fit_l22$Unstand_SE,
    l11_closed_form = l11_cf, l21_closed_form = l21_cf, l22_closed_form = l22_cf,
    share_specific = 1 - rg^2,
    h2_specific = h2_spec, h2_specific_se = h2_spec_se, h2_specific_z = h2_spec / h2_spec_se,
    b = b, neff_subtracted = n1, neff_residual = n2,
    R_effective_fraction = (1 - rg^2) / (1 + b^2 * n2 / n1),
    min_q2_alpha_5e8  = ncp(5e-8)   * sigma_u2,
    min_q2_alpha_25e9 = ncp(2.5e-8) * sigma_u2,
    min_q2_share_of_h2_specific = ncp(2.5e-8) * sigma_u2 / h2_spec
  )
}

chol_summary <- imap_dfr(orderings, summarise_ordering)
glimpse(chol_summary)
write_tsv(ldsc_summary, "ldsc_summary.tsv")
write_tsv(chol_summary, "cholesky_summary.tsv")

## ---- 4. SNP files for userGWAS ---------------------------------------------
## Betas are log-odds from logistic regression -> se.logit = TRUE.
snps <- sumstats(files = unname(clean_files), ref = ref_1000g, trait.names = traits,
                 se.logit = c(TRUE, TRUE), OLS = c(FALSE, FALSE), linprob = c(FALSE, FALSE),
                 maf.filter = maf_filter, keep.indel = FALSE,
                 parallel = n_cores > 1, cores = n_cores)
saveRDS(snps, "sumstats_for_userGWAS.rds")
message("SNPs carried to userGWAS: ", nrow(snps))

writeLines(capture.output(sessionInfo()), "sessionInfo_01.txt")
