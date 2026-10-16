## 04: Miami plot. F1 (shared factor) on top, F2 (sex-specific factor) mirrored below.
suppressPackageStartupMessages({ library(tidyverse) })
source(Sys.getenv("HF_GSEM_CONFIG", "00_config.R"))
setwd(out_dir)
dir.create("plots", showWarnings = FALSE)

y_cap <- NA   # e.g. 30: cap -log10(p) at this value (capped SNPs drawn as diamonds); NA = no cap
## Peaks: triangle on the lead SNP of each peak with p < sig_alpha (hollow = suggestive,
## filled = genome-wide, p < 5e-8). Up for F1, down for F2.
peak_window <- 5e5    # one triangle per peak: the smallest-p SNP, then nothing else within +/- this (bp)

read_factor <- function(label, factor) {
  read_tsv(str_glue("{label}_{factor}_genomewide.tsv.gz"),
           col_select = c(SNP, CHR, BP, Pval_Estimate, pass), show_col_types = FALSE) %>%
    filter(pass, !is.na(Pval_Estimate)) %>%   # pass: no error, S not smoothed (set in 03)
    transmute(SNP, CHR, BP, factor = factor,
              logp = -log10(pmax(Pval_Estimate, .Machine$double.xmin)))
}

## Thin the null bulk without changing the picture: for p > 0.01, keep one point
## per (chromosome, 100 kb bin, 0.01 step in -log10 p).
thin <- function(d) {
  bind_rows(filter(d, logp >= 2),
            filter(d, logp < 2) %>%
              distinct(CHR, bin = BP %/% 1e5, lp = round(logp, 2), .keep_all = TRUE) %>%
              select(-bin, -lp))
}

## Peak markers (display only, distance-based; not an LD clump): take the smallest-p
## SNP, drop SNPs on the same chromosome within +/- peak_window, repeat.
find_peaks <- function(d) {
  rest <- d %>% filter(logp > -log10(sig_alpha)) %>% arrange(desc(logp))
  pk <- rest[0, ]
  while (nrow(rest) > 0) {
    top  <- rest[1, ]
    pk   <- bind_rows(pk, top)
    rest <- filter(rest, !(CHR == top$CHR & abs(BP - top$BP) <= peak_window))
  }
  pk
}

iwalk(orderings, function(ord, label) {
  d <- map_dfr(c("F1", "F2"), ~ read_factor(label, .x))
  n_cap <- 0
  if (!is.na(y_cap)) { n_cap <- sum(d$logp > y_cap); d <- mutate(d, capped = logp > y_cap, logp = pmin(logp, y_cap)) }
  else d <- mutate(d, capped = FALSE)

  ## cumulative genome position
  chr_len <- d %>% group_by(CHR) %>% summarise(len = max(BP), .groups = "drop") %>%
    arrange(CHR) %>% mutate(offset = lag(cumsum(as.numeric(len)), default = 0),
                            mid = offset + len / 2)
  d <- thin(d) %>% left_join(chr_len, by = "CHR") %>%
    mutate(pos = BP + offset,
           y = if_else(factor == "F1", logp, -logp),
           shade = factor(CHR %% 2))

  ymax <- max(abs(d$y), -log10(5e-8)) * 1.15   # headroom for the panel labels
  pk <- d %>% group_by(factor) %>% group_modify(~ find_peaks(.x)) %>% ungroup() %>%
    mutate(y_mark = y + if_else(factor == "F1", 1, -1) * 0.03 * ymax,
           gws = logp > -log10(5e-8))
  write_tsv(pk %>% select(factor, SNP, CHR, BP, logp, genome_wide = gws) %>% mutate(p = 10^-logp),
            str_glue("plots/{label}_miami_peaks.tsv"))
  thr <- tibble(y = c(-log10(5e-8), -log10(sig_alpha)), type = c("5e-8", format(sig_alpha)))
  lab_top <- str_glue("F1: shared factor ({ord[['t1']]} scaled)")
  lab_bot <- str_glue("F2: {ord[['t2']]}-specific factor")

  p <- ggplot(d, aes(pos, y)) +
    geom_point(aes(colour = interaction(factor, shade), shape = capped), size = 0.6, show.legend = FALSE) +
    geom_hline(data = thr, aes(yintercept = y, linetype = type), colour = "grey30", linewidth = 0.3) +
    geom_hline(data = thr, aes(yintercept = -y, linetype = type), colour = "grey30", linewidth = 0.3) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.3) +
    geom_point(data = filter(pk, factor == "F1"), aes(pos, y_mark, fill = gws), shape = 24, size = 1.8, show.legend = FALSE) +
    geom_point(data = filter(pk, factor == "F2"), aes(pos, y_mark, fill = gws), shape = 25, size = 1.8, show.legend = FALSE) +
    scale_fill_manual(values = c(`TRUE` = "black", `FALSE` = "white")) +
    annotate("text", x = 0, y =  ymax, label = lab_top, hjust = 0, vjust = 1, size = 3.5) +
    annotate("text", x = 0, y = -ymax, label = lab_bot, hjust = 0, vjust = 0, size = 3.5) +
    scale_colour_manual(values = c(F1.0 = "#1f4e79", F1.1 = "#6fa8dc",
                                   F2.0 = "#7f2704", F2.1 = "#f4a261")) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18)) +
    scale_linetype_manual(values = setNames(c("dashed", "dotted"), thr$type), name = "p threshold") +
    scale_x_continuous(breaks = chr_len$mid, labels = chr_len$CHR, expand = expansion(mult = 0.005)) +
    scale_y_continuous(labels = abs, limits = c(-ymax, ymax)) +
    labs(x = "Chromosome", y = expression(-log[10](italic(p))),
         caption = if (n_cap > 0) str_glue("{n_cap} SNP rows capped at -log10(p) = {y_cap} (diamonds)") else NULL) +
    theme_classic(base_size = 11) +
    theme(legend.position = "top", panel.grid = element_blank())

  ggsave(str_glue("plots/{label}_miami.png"), p, width = 12, height = 6, dpi = 300)
  ggsave(str_glue("plots/{label}_miami.pdf"), p, width = 12, height = 6)
  message(str_glue("{label}: plotted {nrow(d)} points (after thinning); peaks F1 = {sum(pk$factor == 'F1')}, F2 = {sum(pk$factor == 'F2')}"))
})
