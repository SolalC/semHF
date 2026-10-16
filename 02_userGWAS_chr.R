## 02: userGWAS for one chromosome (SLURM array task = chromosome).
suppressPackageStartupMessages({ library(GenomicSEM); library(tidyverse) })
source(Sys.getenv("HF_GSEM_CONFIG", "00_config.R"))
set.seed(seed)
setwd(out_dir)

chr <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", NA))
if (is.na(chr)) chr <- as.integer(commandArgs(trailingOnly = TRUE)[1])
stopifnot(chr %in% 1:22)

covstruc <- readRDS("covstruc.rds")
snps     <- readRDS("sumstats_for_userGWAS.rds") %>% filter(CHR == chr)
message("Chromosome ", chr, ": ", nrow(snps), " SNPs; cores = ", n_cores)
dir.create("userGWAS", showWarnings = FALSE)

iwalk(orderings, function(ord, label) {
  res <- userGWAS(covstruc = covstruc, SNPs = snps, estimation = "DWLS",
                  model = cholesky_model(ord[["t1"]], ord[["t2"]], with_snp = TRUE),
                  sub = c("F1~SNP", "F2~SNP"),
                  parallel = n_cores > 1, cores = n_cores,
                  fix_measurement = TRUE, smooth_check = TRUE)
  ## res[[1]] = F1 ~ SNP (shared factor), res[[2]] = F2 ~ SNP (sex-specific factor)
  ## lavaan warnings contain newlines/tabs, which write_tsv does not quote: flatten them.
  res <- map(res, ~ mutate(.x, across(where(is.character), \(v) str_squish(str_replace_all(v, "[\t\r\n]+", " ")))))
  write_tsv(res[[1]], file.path("userGWAS", str_glue("{label}_F1_chr{chr}.tsv.gz")))
  write_tsv(res[[2]], file.path("userGWAS", str_glue("{label}_F2_chr{chr}.tsv.gz")))
})

writeLines(capture.output(sessionInfo()), str_glue("userGWAS/sessionInfo_02_chr{chr}.txt"))
