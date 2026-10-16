## 03: combine chromosomes, QC, and check userGWAS against the closed form
##     z_F2 = (beta_2 - b beta_1) / sqrt(se_2^2 + b^2 se_1^2)   (no sample overlap)
suppressPackageStartupMessages({ library(tidyverse) })
source(Sys.getenv("HF_GSEM_CONFIG", "00_config.R"))
setwd(out_dir)

chol_summary <- read_tsv("cholesky_summary.tsv", show_col_types = FALSE)
snps <- readRDS("sumstats_for_userGWAS.rds")

## userGWAS writes lavaan warnings that contain newlines into the `warning`/`error`
## columns, and write_tsv does not quote them, so one SNP can span several lines.
## Repair: a record starts on a line whose first two fields are a known SNP and its
## chromosome; any other line is a continuation and is glued back onto its record.
read_userGWAS_file <- function(path, chr, snp_ids) {
  lines  <- str_remove_all(read_lines(path), "\r")
  header <- lines[1]; body <- lines[-1]
  f1 <- str_extract(body, "^[^\t]*")
  f2 <- str_match(body, "^[^\t]*\t([^\t]*)")[, 2]
  start <- f1 %in% snp_ids & f2 == as.character(chr)
  starts <- which(start); conts <- which(!start)
  conts_orphan <- conts[conts < min(c(starts, Inf))]
  conts <- setdiff(conts, conts_orphan)
  if (length(conts) > 0) {
    owner <- starts[findInterval(conts, starts)]
    glued <- tapply(body[conts], owner, str_c, collapse = " ")
    i <- as.integer(names(glued))
    body[i] <- str_c(body[i], " ", glued)
  }
  d <- read_tsv(I(c(header, body[starts])), col_types = cols(.default = col_character()))
  attr(d, "repair") <- tibble(file = basename(path), n_records = length(starts),
                              n_continuation_lines = length(conts), n_orphan_lines = length(conts_orphan),
                              n_snps_expected = sum(snps$CHR == chr))
  d
}

## Read every column as character, bind, then convert types once: per-file type
## guessing differs between chromosomes (e.g. CHR read as character in one file).
combine_factor <- function(label, factor) {
  files <- list.files("userGWAS", pattern = str_glue("^{label}_{factor}_chr[0-9]+\\.tsv\\.gz$"), full.names = TRUE)
  chrs  <- as.integer(str_match(basename(files), "_chr([0-9]+)\\.tsv")[, 2])
  d <- map2(files, chrs, read_userGWAS_file, snp_ids = snps$SNP)
  repair <- map_dfr(d, attr, "repair")
  write_tsv(repair, str_glue("{label}_{factor}_rows_per_chr.tsv"))
  print(repair)
  if (length(files) != 22 || any(repair$n_records != repair$n_snps_expected)) {
    warning(str_glue("{label} {factor}: {length(files)} chromosome files; record counts differ from ",
                     "the SNP counts for some chromosomes. See {label}_{factor}_rows_per_chr.tsv."))
  }
  d <- bind_rows(d)

  ## Force numeric result columns; record values that are not numbers (e.g. text
  ## written by userGWAS for SNPs that failed) instead of letting them break the QC.
  num_cols <- intersect(c("MAF", "est", "se_c", "SE", "Z_Estimate", "Pval_Estimate",
                          "chisq", "chisq_df", "chisq_pval", "AIC", "Z_smooth"), names(d))
  bad <- map_dfr(num_cols, function(col) {
    x <- d[[col]]
    is_bad <- !is.na(x) & x != "NA" & is.na(suppressWarnings(as.numeric(x)))
    tibble(column = col, n_non_numeric = sum(is_bad), examples = str_c(head(unique(x[is_bad]), 5), collapse = " | "))
  })
  write_tsv(bad, str_glue("{label}_{factor}_non_numeric_values.tsv"))
  if (any(bad$n_non_numeric > 0)) print(bad)

  d %>%
    mutate(across(all_of(num_cols), ~ suppressWarnings(as.numeric(.x))),
           CHR = as.integer(CHR), BP = as.numeric(BP)) %>%
    ## pass = no lavaan error and S not smoothed. A smoothed S means the SNP's input
    ## effect was too large for the trait h2 (S + SNP not positive definite); the fit
    ## is then to an altered matrix and its est/SE are not interpretable.
    mutate(pass = (is.na(error) | error == "0") &
             (if ("Z_smooth" %in% names(.)) is.na(Z_smooth) | Z_smooth <= z_smooth_max else TRUE)) %>%
    arrange(CHR, BP)
}

## SNPs below sig_alpha (no clumping; do LD clumping separately with a reference panel).
extract_significant <- function(d, label, factor, alpha = sig_alpha) {
  sig <- d %>%
    filter(Pval_Estimate < alpha, pass) %>%
    left_join(snps %>% select(SNP, starts_with("beta."), starts_with("se.")), by = "SNP") %>%
    arrange(Pval_Estimate)
  tag <- format(alpha, scientific = TRUE)
  write_tsv(sig, str_glue("{label}_{factor}_significant_p{tag}.tsv"))
  nrow(sig)
}

qc <- imap_dfr(orderings, function(ord, label) {
  t1 <- ord[["t1"]]; t2 <- ord[["t2"]]
  pars <- chol_summary %>% filter(ordering == label)

  out <- map(c(F1 = "F1", F2 = "F2"), ~ combine_factor(label, .x))
  iwalk(out, ~ write_tsv(.x, str_glue("{label}_{.y}_genomewide.tsv.gz")))
  n_sig <- imap(out, ~ extract_significant(.x, label, .y))

  ## Closed-form z for the sex-specific factor, from the standardised betas
  cf <- snps %>%
    transmute(SNP,
              z_closed_form = (.data[[str_c("beta.", t2)]] - pars$b * .data[[str_c("beta.", t1)]]) /
                sqrt(.data[[str_c("se.", t2)]]^2 + pars$b^2 * .data[[str_c("se.", t1)]]^2))
  chk <- out$F2 %>% filter(pass) %>% select(SNP, Z_Estimate) %>% inner_join(cf, by = "SNP") %>%
    filter(is.finite(Z_Estimate), is.finite(z_closed_form))

  imap_dfr(out, function(d_all, f) {
    d <- filter(d_all, pass)
    tibble(ordering = label, factor = f, n_snps = nrow(d_all),
           n_error = if ("error" %in% names(d_all)) sum(d_all$error != 0, na.rm = TRUE) else NA_integer_,
           n_smoothed_dropped = if ("Z_smooth" %in% names(d_all)) sum(d_all$Z_smooth > z_smooth_max, na.rm = TRUE) else NA_integer_,
           n_pass = nrow(d),
           lambda_gc = median(d$Z_Estimate^2, na.rm = TRUE) / qchisq(0.5, 1),
           n_p_5e8 = sum(d$Pval_Estimate < 5e-8, na.rm = TRUE),
           sig_alpha = sig_alpha,
           n_p_sig_alpha = n_sig[[f]])
  }) %>%
    mutate(cor_z_userGWAS_vs_closed_form = if_else(factor == "F2", cor(chk$Z_Estimate, chk$z_closed_form), NA_real_))
})

print(qc)
write_tsv(qc, "userGWAS_qc.tsv")
writeLines(capture.output(sessionInfo()), "sessionInfo_03.txt")
