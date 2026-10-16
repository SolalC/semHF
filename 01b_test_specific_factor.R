## 01b (optional): chi-square difference test of the female-specific path.
## Wald z for l22 (a square root near its boundary) is not reliable; this fits the
## Cholesky model with l22 fixed to 0 (df = 1) against the saturated model (df = 0).
suppressPackageStartupMessages({ library(GenomicSEM); library(tidyverse) })
source(Sys.getenv("HF_GSEM_CONFIG", "00_config.R"))
setwd(out_dir)

covstruc <- readRDS("covstruc.rds")
ord <- orderings[["F_specific"]]
## Setting l22 = 0 removes F2 entirely: drop every F2 line from the Cholesky model.
## The result is one factor with both residual variances fixed to 0, i.e. r_g = 1.
model_l22_zero <- str_split(cholesky_model(ord[["t1"]], ord[["t2"]]), "\n")[[1]] %>%
  str_trim() %>%
  discard(~ str_detect(.x, "F2")) %>%
  str_c(collapse = "\n")
cat(model_l22_zero, "\n")

fit0 <- usermodel(covstruc, estimation = "DWLS", model = model_l22_zero, CFIcalc = FALSE)
print(fit0$modelfit)
## Saturated model has chisq = 0, df = 0, so the difference test is fit0's chisq on 1 df.
