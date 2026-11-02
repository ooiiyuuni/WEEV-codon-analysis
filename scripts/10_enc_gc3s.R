# ============================================================
# 10_enc_gc3s.R
# Purpose : ENc-GC3s plots (full range and zoom) for WEEV CDS, NSP and SP,
#           colored by host group; compute expected ENc (ENCe), ENC ratio
#           (ENCr) and a one-sample Wilcoxon signed-rank test of ENCr vs 0.
# Input   : results/tables/codonw/CodonUsageIndices_CDS.csv
#           results/tables/codonw/CodonUsageIndices_NSP.csv
#           results/tables/codonw/CodonUsageIndices_SP.csv
#           results/tables/metadata.csv
# Output  : results/figures/ENc_GC3s_all.pdf
#           results/tables/ENc_GC3s/ENc_GC3s_per_sequence_CDS.csv
#           results/tables/ENc_GC3s/ENc_GC3s_per_sequence_NSP.csv
#           results/tables/ENc_GC3s/ENc_GC3s_per_sequence_SP.csv
#           results/tables/ENc_GC3s/ENc_GC3s_wilcox.csv
# Dependencies: here, tidyverse, patchwork
# Random seed : none (no stochastic step; Wilcoxon with exact = FALSE is deterministic)
# Run     : Rscript scripts/10_enc_gc3s.R
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(patchwork)
})

fam <- "ArialMT"

tab_dir <- here("results", "tables")
fig_dir <- here("results", "figures")
out_tab <- here("results", "tables", "ENc_GC3s")
dir.create(out_tab, showWarnings = FALSE, recursive = TRUE)

meta <- read.csv(file.path(tab_dir, "metadata.csv"),
                 check.names = FALSE, stringsAsFactors = FALSE) %>%
  select(Accession, Host) %>%
  mutate(meta_row = row_number())

host_group <- function(h) {
  h <- tolower(trimws(h))
  dplyr::case_when(
    h %in% c("equine", "horse") ~ "Equine",
    h %in% c("culex tarsalis", "mosquito", "culex sp.") ~ "Mosquito",
    h %in% c("homo sapiens") ~ "Human",
    TRUE ~ "Others"
  )
}

read_region <- function(region) {
  ind <- read.csv(file.path(tab_dir, "codonw", paste0("CodonUsageIndices_", region, ".csv")),
                  check.names = FALSE, stringsAsFactors = FALSE) %>%
    filter(!grepl("_Mean$|_SD$", Accession)) %>%
    select(Accession, ENCo = Nc, GC3s)
  ind %>%
    inner_join(meta, by = "Accession") %>%
    mutate(Region = region,
           Group  = factor(host_group(Host),
                           levels = c("Mosquito", "Equine", "Human", "Others"))) %>%
    mutate(ENCe = 2 + GC3s + 29 / (GC3s^2 + (1 - GC3s)^2),
           ENCr = (ENCe - ENCo) / ENCe)
}

regions <- c("CDS", "NSP", "SP")
dat_all <- bind_rows(lapply(regions, read_region))

# ---- sanity checks ----
stopifnot(nrow(dat_all) == 156)
for (r in regions) {
  d <- dat_all %>% filter(Region == r)
  stopifnot(nrow(d) == 52)
  stopifnot(!any(is.na(d$Host)), !any(is.na(d$GC3s)), !any(is.na(d$ENCo)))
}
ord <- lapply(regions, function(r) dat_all %>% filter(Region == r) %>% pull(Accession))
stopifnot(identical(ord[[1]], ord[[2]]), identical(ord[[1]], ord[[3]]))

# ---- theoretical curve ----
make_curve <- function(xr) {
  data.frame(GC3s = seq(xr[1], xr[2], length.out = 500)) %>%
    mutate(ENCe = 2 + GC3s + 29 / (GC3s^2 + (1 - GC3s)^2))
}

x_all <- range(dat_all$GC3s)
y_all <- range(dat_all$ENCo)

x_full <- c(0.45, 0.55)
curve_full <- make_curve(x_full)
y_full <- range(c(dat_all$ENCo, curve_full$ENCe))
y_full <- y_full + c(-1, 1) * diff(y_full) * 0.05

x_zoom <- c(0.485, x_all[2] + diff(x_all) * 0.05)
y_zoom <- y_all + c(-1, 1) * diff(y_all) * 0.05

grp_cols <- c(Mosquito = "#64A9E2", Equine = "#DDCE43",
              Human = "#E41A1C", Others = "#681F72")

base_theme <- theme_classic(base_size = 10, base_family = fam) +
  theme(
    plot.title   = element_text(size = 16, face = "bold", hjust = 0.5,
                                family = fam, margin = margin(b = 8)),
    axis.title.x = element_text(size = 16, face = "bold", family = fam,
                                margin = margin(t = 10)),
    axis.title.y = element_text(size = 16, face = "bold", family = fam,
                                margin = margin(r = 20)),
    axis.text    = element_text(size = 10, color = "black", family = fam),
    axis.line    = element_line(color = "black", linewidth = 0.75 * 0.352778),
    axis.ticks   = element_line(color = "black", linewidth = 0.75 * 0.352778),
    legend.title = element_text(size = 10, face = "bold", family = fam),
    legend.text  = element_text(size = 10, family = fam),
    plot.margin  = margin(8, 8, 8, 8),
    panel.grid.minor = element_blank()
  )

grid_theme <- theme(
  panel.grid.major = element_line(color = "grey85", linewidth = 0.3)
)

make_plot_full <- function(df, curve_df, x_rng, y_rng, title, show_legend = FALSE) {
  p <- ggplot() +
    geom_line(data = curve_df, aes(x = GC3s, y = ENCe),
              color = "black", linewidth = 0.8) +
    geom_point(data = df, aes(x = GC3s, y = ENCo, color = Group),
               shape = 16, size = 2.4, alpha = 0.9) +
    scale_color_manual(values = grp_cols, name = "Host group") +
    scale_x_continuous(limits = x_rng, expand = c(0, 0),
                       breaks = seq(x_rng[1], x_rng[2], by = 0.05)) +
    scale_y_continuous(limits = y_rng, expand = c(0, 0),
                       breaks = scales::breaks_pretty(n = 5)) +
    coord_fixed(ratio = (x_rng[2] - x_rng[1]) / (y_rng[2] - y_rng[1])) +
    labs(title = title, x = "GC3s", y = "Effective number of codons (ENc)") +
    base_theme + theme(panel.grid = element_blank())
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

make_plot_zoom <- function(df, x_rng, y_rng, title, show_legend = FALSE) {
  p <- ggplot() +
    geom_point(data = df, aes(x = GC3s, y = ENCo, color = Group),
               shape = 16, size = 2.4, alpha = 0.9) +
    scale_color_manual(values = grp_cols, name = "Host group") +
    scale_x_continuous(limits = x_rng, expand = c(0, 0),
                       breaks = seq(0.49, 0.53, by = 0.01)) +
    scale_y_continuous(limits = y_rng, expand = c(0, 0),
                       breaks = scales::breaks_pretty(n = 5)) +
    coord_fixed(ratio = (x_rng[2] - x_rng[1]) / (y_rng[2] - y_rng[1])) +
    labs(title = title, x = "GC3s", y = "Effective number of codons (ENc)") +
    base_theme + grid_theme
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

plots <- list()
for (r in regions) {
  df_r <- dat_all %>% filter(Region == r)
  plots[[paste0(r, "_full")]] <- make_plot_full(df_r, curve_full, x_full, y_full,
                                                paste0(r, " (full range)"),
                                                show_legend = (r == "CDS"))
  plots[[paste0(r, "_zoom")]] <- make_plot_zoom(df_r, x_zoom, y_zoom,
                                                paste0(r, " (zoom)"),
                                                show_legend = FALSE)
}

combined <- (plots[["CDS_full"]] | plots[["NSP_full"]] | plots[["SP_full"]]) /
            (plots[["CDS_zoom"]] | plots[["NSP_zoom"]] | plots[["SP_zoom"]]) +
  plot_layout(widths = c(1, 1, 1), heights = c(1, 1), guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "ENc_GC3s_all.pdf"),
       combined, width = 18, height = 10, device = "pdf")

# ---- per-sequence tables (one per region) ----
for (r in regions) {
  df_r <- dat_all %>% filter(Region == r)
  seq_out <- df_r %>%
    select(Accession, Host, Group, GC3s, ENCo, ENCe, ENCr, meta_row) %>%
    arrange(meta_row) %>%
    select(-meta_row)
  write.csv(seq_out, file.path(out_tab, paste0("ENc_GC3s_per_sequence_", r, ".csv")), row.names = FALSE)
}

# ---- single combined Wilcoxon table (3 rows) ----
wilcox_all <- bind_rows(lapply(regions, function(r) {
  df_r <- dat_all %>% filter(Region == r)
  wt <- wilcox.test(df_r$ENCr, mu = 0, alternative = "two.sided", exact = FALSE)
  data.frame(
    region = r, n = nrow(df_r),
    median_ENCr = median(df_r$ENCr), mean_ENCr = mean(df_r$ENCr),
    n_ENCr_gt_0.05 = sum(df_r$ENCr > 0.05),
    prop_ENCr_gt_0.05 = mean(df_r$ENCr > 0.05),
    statistic = unname(wt$statistic), p_value = wt$p.value
  )
}))
write.csv(wilcox_all, file.path(out_tab, "ENc_GC3s_wilcox.csv"), row.names = FALSE)

cat("Done.\n")
cat("Figure:", file.path(fig_dir, "ENc_GC3s_all.pdf"), "\n")
cat("Tables:", out_tab, "\n")
