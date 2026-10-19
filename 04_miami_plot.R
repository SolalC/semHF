## 04: Miami plots, all on the same y-scale.
##   {ordering}_miami : F1 (shared factor) on top, F2 (sex-specific factor) below, per ordering
##   F1_miami, F2_miami : the same factor, female-specific ordering on top, male-specific below
suppressPackageStartupMessages({ library(tidyverse) })
source(Sys.getenv("HF_GSEM_CONFIG", "00_config.R"))
setwd(out_dir)
dir.create("plots", showWarnings = FALSE)

y_cap <- NA   # e.g. 30: cap -log10(p) at this value (capped SNPs drawn as diamonds); NA = no cap
## Peaks: triangle on the lead SNP of each peak with p < sig_alpha (hollow = suggestive,
## filled = genome-wide, p < 5e-8). Up for the top panel, down for the bottom panel.
peak_window <- 5e5    # one triangle per peak: the smallest-p SNP, then nothing else within +/- this (bp)

read_factor <- function(label, factor) {
  read_tsv(str_glue("{label}_{factor}_genomewide.tsv.gz"),
           col_select = c(SNP, CHR, BP, Pval_Estimate, pass), show_col_types = FALSE) %>%
    filter(pass, !is.na(Pval_Estimate)) %>%   # pass: no error, S not smoothed (set in 03)
    transmute(SNP, CHR, BP, logp = -log10(pmax(Pval_Estimate, .Machine$double.xmin)))
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

## ---- read everything once ---------------------------------------------------
res <- expand_grid(ordering = names(orderings), factor = c("F1", "F2")) %>%
  mutate(data = map2(ordering, factor, read_factor))
all_d <- bind_rows(res$data)
if (!is.na(y_cap)) message(sum(all_d$logp > y_cap), " SNP rows above y_cap = ", y_cap)

## Shared x-axis (cumulative position) and shared y-scale for every plot
chr_len <- all_d %>% group_by(CHR) %>% summarise(len = max(BP), .groups = "drop") %>%
  arrange(CHR) %>% mutate(offset = lag(cumsum(as.numeric(len)), default = 0),
                          mid = offset + len / 2)
ymax <- max(if (is.na(y_cap)) max(all_d$logp) else min(max(all_d$logp), y_cap), -log10(5e-8)) * 1.15
rm(all_d)

get_d <- function(ord, f) res$data[[which(res$ordering == ord & res$factor == f)]]

panel_label <- function(ord, f) {
  o <- orderings[[ord]]
  if (f == "F1") str_glue("{ord} F1: shared factor ({o[['t1']]} scaled)")
  else           str_glue("{ord} F2: {o[['t2']]}-specific factor")
}

## ---- one Miami plot -----------------------------------------------------------
miami <- function(top, bottom, lab_top, lab_bot, file) {
  d <- bind_rows(top = top, bottom = bottom, .id = "panel")
  if (!is.na(y_cap)) d <- mutate(d, capped = logp > y_cap, logp = pmin(logp, y_cap))
  else d <- mutate(d, capped = FALSE)
  n_cap <- sum(d$capped)

  d <- d %>% group_by(panel) %>% group_modify(~ thin(.x)) %>% ungroup() %>%
    left_join(chr_len, by = "CHR") %>%
    mutate(pos = BP + offset,
           y = if_else(panel == "top", logp, -logp),
           shade = factor(CHR %% 2))

  pk <- d %>% group_by(panel) %>% group_modify(~ find_peaks(.x)) %>% ungroup() %>%
    mutate(y_mark = y + if_else(panel == "top", 1, -1) * 0.03 * ymax,
           gws = logp > -log10(5e-8))
  write_tsv(pk %>% transmute(panel = if_else(panel == "top", lab_top, lab_bot),
                             SNP, CHR, BP, logp, p = 10^-logp, genome_wide = gws),
            str_glue("plots/{file}_peaks.tsv"))
  thr <- tibble(y = c(-log10(5e-8), -log10(sig_alpha)), type = c("5e-8", format(sig_alpha)))

  p <- ggplot(d, aes(pos, y)) +
    geom_point(aes(colour = interaction(panel, shade), shape = capped), size = 0.6, show.legend = FALSE) +
    geom_hline(data = thr, aes(yintercept = y, linetype = type), colour = "grey30", linewidth = 0.3) +
    geom_hline(data = thr, aes(yintercept = -y, linetype = type), colour = "grey30", linewidth = 0.3) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.3) +
    geom_point(data = filter(pk, panel == "top"),    aes(pos, y_mark, fill = gws), shape = 24, size = 1.8, show.legend = FALSE) +
    geom_point(data = filter(pk, panel == "bottom"), aes(pos, y_mark, fill = gws), shape = 25, size = 1.8, show.legend = FALSE) +
    scale_fill_manual(values = c(`TRUE` = "black", `FALSE` = "white")) +
    annotate("text", x = 0, y =  ymax, label = lab_top, hjust = 0, vjust = 1, size = 3.5) +
    annotate("text", x = 0, y = -ymax, label = lab_bot, hjust = 0, vjust = 0, size = 3.5) +
    scale_colour_manual(values = c(top.0 = "#1f4e79", top.1 = "#6fa8dc",
                                   bottom.0 = "#7f2704", bottom.1 = "#f4a261")) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18)) +
    scale_linetype_manual(values = setNames(c("dashed", "dotted"), thr$type), name = "p threshold") +
    scale_x_continuous(breaks = chr_len$mid, labels = chr_len$CHR, expand = expansion(mult = 0.005)) +
    scale_y_continuous(labels = abs, limits = c(-ymax, ymax)) +
    labs(x = "Chromosome", y = expression(-log[10](italic(p))),
         caption = if (n_cap > 0) str_glue("{n_cap} SNP rows capped at -log10(p) = {y_cap} (diamonds)") else NULL) +
    theme_classic(base_size = 11) +
    theme(legend.position = "top", panel.grid = element_blank())

  ggsave(str_glue("plots/{file}.png"), p, width = 12, height = 6, dpi = 300)
  ggsave(str_glue("plots/{file}.pdf"), p, width = 12, height = 6)
  message(str_glue("{file}: {nrow(d)} points after thinning; peaks top = {sum(pk$panel == 'top')}, ",
                   "bottom = {sum(pk$panel == 'bottom')}"))
}

## ---- plots ----------------------------------------------------------------------
## Per ordering: F1 on top, F2 below
for (ord in names(orderings)) {
  miami(get_d(ord, "F1"), get_d(ord, "F2"), panel_label(ord, "F1"), panel_label(ord, "F2"),
        str_glue("{ord}_miami"))
}

## Per factor: female-specific ordering on top, male-specific below
if (all(c("F_specific", "M_specific") %in% names(orderings))) {
  for (f in c("F1", "F2")) {
    miami(get_d("F_specific", f), get_d("M_specific", f),
          panel_label("F_specific", f), panel_label("M_specific", f), str_glue("{f}_miami"))
  }
}
message("All plots share the y-axis: 0 to ", round(ymax, 1))
