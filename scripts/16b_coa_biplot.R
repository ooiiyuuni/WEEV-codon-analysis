# ============================================================
# 16b_coa_biplot.R
# Purpose : COA biplot for CDS, NSP and SP, colored by host group
#           and by geographic region. Axes are drawn as arrows,
#           tick marks are removed, grey major grid is kept, and
#           all six panels share the same symmetric limits with
#           the origin at the centre. Variable arrows are drawn
#           from the Spearman correlation of each index with
#           Axis 1 and Axis 2.
# Input   : results/tables/COA/COA_coords_{CDS,NSP,SP}.csv
#           results/tables/COA/COA_axis_contribution.csv
#           results/tables/codonw/CodonUsageIndices_{CDS,NSP,SP}.csv
#           results/tables/BaseComposition_{CDS,NSP,SP}.csv
#           results/tables/metadata.csv
# Output  : results/figures/COA_biplot_by_host_{CDS,NSP,SP}.pdf
#           results/figures/COA_biplot_by_region_{CDS,NSP,SP}.pdf
#           results/tables/COA/COA_arrow_correlations.csv
#           results/tables/COA/COA_axis_contribution_wide.csv
# Dependencies: here, tidyverse
# Random seed : not applicable (no stochastic step)
# Run     : Rscript scripts/16b_coa_biplot.R
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
})

coa_dir <- here("results", "tables", "COA")
out_fig <- here("results", "figures")
if (!dir.exists(out_fig)) dir.create(out_fig, recursive = TRUE)

coords <- map_dfr(c("CDS", "NSP", "SP"), function(r) {
  read_csv(file.path(coa_dir, paste0("COA_coords_", r, ".csv")), show_col_types = FALSE) %>%
    mutate(Region = r)
}) %>%
  mutate(Region = factor(Region, levels = c("CDS", "NSP", "SP")))

meta <- read_csv(here("results", "tables", "metadata.csv"), show_col_types = FALSE) %>%
  select(Accession, Host, Country)

host_group <- function(x) {
  x_low <- tolower(x)
  case_when(
    x_low %in% c("culex tarsalis", "mosquito", "culex sp.") ~ "Mosquito",
    x_low %in% c("equine", "horse") ~ "Equine",
    x_low %in% c("homo sapiens") ~ "Human",
    TRUE ~ "Others"
  )
}

region_group <- function(x) {
  prefix <- str_extract(x, "^[^:]+") %>% str_trim()
  case_when(
    prefix %in% c("Uruguay", "Brazil", "Argentina") ~ "South America",
    prefix %in% c("USA", "Canada") ~ "North America",
    TRUE ~ NA_character_
  )
}

meta <- meta %>%
  mutate(
    HostGroup = factor(host_group(Host),
                       levels = c("Mosquito", "Equine", "Human", "Others")),
    GeoRegion = factor(region_group(Country),
                       levels = c("North America", "South America"))
  )
stopifnot(!any(is.na(meta$HostGroup)))
stopifnot(!any(is.na(meta$GeoRegion)))

read_vars <- function(region) {
  cod <- read_csv(here("results", "tables", "codonw",
                       paste0("CodonUsageIndices_", region, ".csv")),
                  show_col_types = FALSE) %>%
    filter(!str_detect(Accession, "_Mean$|_SD$")) %>%
    select(Accession, A3s, U3s, G3s, C3s, GC3s, Nc, Gravy, Aromo)

  bc <- read_csv(here("results", "tables",
                      paste0("BaseComposition_", region, ".csv")),
                 show_col_types = FALSE) %>%
    filter(!str_detect(Accession, "_Mean$|_SD$")) %>%
    select(Accession, `A%`, `U%`, `G%`, `C%`, `GC%`) %>%
    rename(A = `A%`, U = `U%`, G = `G%`, C = `C%`, GC = `GC%`)

  inner_join(cod, bc, by = "Accession")
}

vars_by_region <- lapply(c(CDS = "CDS", NSP = "NSP", SP = "SP"), read_vars)

for (nm in names(vars_by_region)) {
  stopifnot(nrow(vars_by_region[[nm]]) == 52)
  stopifnot(!any(is.na(vars_by_region[[nm]])))
}

dat <- map_dfr(c("CDS", "NSP", "SP"), function(r) {
  coords %>% filter(Region == r) %>%
    inner_join(vars_by_region[[r]], by = "Accession")
}) %>%
  left_join(meta, by = "Accession") %>%
  mutate(Region = factor(Region, levels = c("CDS", "NSP", "SP")))

stopifnot(nrow(dat) == 52 * 3)
stopifnot(!any(is.na(dat$HostGroup)))
stopifnot(!any(is.na(dat$GeoRegion)))

var_names <- c("GC3s", "Nc", "Gravy", "Aromo",
               "A", "U", "G", "C", "GC",
               "A3s", "U3s", "G3s", "C3s")

arrow_tab <- map_dfr(c("CDS", "NSP", "SP"), function(r) {
  d <- dat %>% filter(Region == r)
  map_dfr(var_names, function(v) {
    x <- d[[v]]
    ct1 <- suppressWarnings(cor.test(x, d$Axis1, method = "spearman"))
    ct2 <- suppressWarnings(cor.test(x, d$Axis2, method = "spearman"))
    tibble(
      Region = r,
      Variable = v,
      r_Axis1 = unname(ct1$estimate),
      P_Axis1 = ct1$p.value,
      r_Axis2 = unname(ct2$estimate),
      P_Axis2 = ct2$p.value
    )
  })
})
write_csv(arrow_tab, file.path(coa_dir, "COA_arrow_correlations.csv"))

contrib <- read_csv(file.path(coa_dir, "COA_axis_contribution.csv"), show_col_types = FALSE)

axis_lab <- function(region) {
  v <- contrib %>% filter(Region == !!region, Axis %in% c(1, 2)) %>% arrange(Axis)
  paste0("Axis ", v$Axis, " (", sprintf("%.1f", v$Contribution_percent), "%)")
}

axis_wide <- contrib %>%
  filter(Axis <= 4) %>%
  select(Region, Axis, Contribution_percent, Cumulative_percent) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  pivot_wider(
    names_from = Region,
    values_from = c(Contribution_percent, Cumulative_percent),
    names_glue = "{Region}_{.value}"
  ) %>%
  select(Axis,
         starts_with("CDS_"), starts_with("NSP_"), starts_with("SP_"))

write_csv(axis_wide, file.path(coa_dir, "COA_axis_contribution_wide.csv"))

r_x <- max(abs(dat$Axis1)) * 1.08
r_y <- max(abs(dat$Axis2)) * 1.08
x_lim <- c(-r_x, r_x)
y_lim <- c(-r_y, r_y)
br_x <- pretty(x_lim, n = 5)
br_y <- pretty(y_lim, n = 5)

arr_x <- r_x * 0.8
arr_y <- r_y * 0.8

host_cols <- c(Mosquito = "#64A9E2", Equine = "#DDCE43",
               Human = "#E41A1C", Others = "#681F72")
geo_cols  <- c(`North America` = "#2CA25F", `South America` = "#F16913")

common_theme <- theme_classic(base_size = 10, base_family = "ArialMT") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0, family = "ArialMT"),
    axis.text = element_text(size = 10, color = "black", family = "ArialMT"),
    axis.title = element_text(size = 16, face = "bold", family = "ArialMT"),
    axis.title.y = element_text(margin = margin(r = 20)),
    axis.title.x = element_text(margin = margin(t = 12)),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 10, family = "ArialMT"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

make_plot <- function(region, colour_by) {
  d <- dat %>% filter(Region == region)
  a <- arrow_tab %>%
    filter(Region == !!region) %>%
    mutate(xend = r_Axis1 * arr_x, yend = r_Axis2 * arr_y)
  xl <- axis_lab(region)[1]
  yl <- axis_lab(region)[2]

  if (colour_by == "host") {
    p <- ggplot(d, aes(Axis1, Axis2, colour = HostGroup)) +
      geom_point(size = 2.2, alpha = 0.9) +
      scale_colour_manual(values = host_cols, drop = FALSE)
  } else {
    p <- ggplot(d, aes(Axis1, Axis2, colour = GeoRegion)) +
      geom_point(size = 2.2, alpha = 0.9) +
      scale_colour_manual(values = geo_cols, drop = FALSE)
  }

  p +
    annotate("segment",
             x = -r_x, y = 0, xend = r_x, yend = 0,
             arrow = arrow(length = unit(0.02, "npc")),
             colour = "black", linewidth = 0.75 * 0.352778) +
    annotate("segment",
             x = 0, y = -r_y, xend = 0, yend = r_y,
             arrow = arrow(length = unit(0.02, "npc")),
             colour = "black", linewidth = 0.75 * 0.352778) +
    geom_segment(data = a,
                 aes(x = 0, y = 0, xend = xend, yend = yend),
                 inherit.aes = FALSE,
                 arrow = arrow(length = unit(0.02, "npc")),
                 colour = "grey20", linewidth = 0.4) +
    geom_text(data = a,
              aes(x = xend * 1.06, y = yend * 1.06, label = Variable),
              inherit.aes = FALSE,
              size = 10 * 0.352778, family = "ArialMT", colour = "grey20") +
    scale_x_continuous(limits = x_lim, breaks = br_x, expand = c(0, 0)) +
    scale_y_continuous(limits = y_lim, breaks = br_y, expand = c(0, 0)) +
    labs(title = paste0(region, " - ",
                        ifelse(colour_by == "host", "by host group", "by geographic region")),
         x = xl, y = yl) +
    common_theme
}

# ---------- combine and export figures ----------
library(patchwork)

# 1. Host plots
p_host_CDS <- make_plot("CDS", "host")
p_host_NSP <- make_plot("NSP", "host")
p_host_SP  <- make_plot("SP", "host")
combined_host <- p_host_CDS + p_host_NSP + p_host_SP +
  plot_layout(guides = "collect") & theme(legend.position = "right")
ggsave(file.path(out_fig, "COA_biplot_by_host_all.pdf"),
       combined_host, width = 18, height = 5, device = "pdf")

# 2. Geographic region plots
p_geo_CDS <- make_plot("CDS", "region")
p_geo_NSP <- make_plot("NSP", "region")
p_geo_SP  <- make_plot("SP", "region")
combined_geo <- p_geo_CDS + p_geo_NSP + p_geo_SP +
  plot_layout(guides = "collect") & theme(legend.position = "right")
ggsave(file.path(out_fig, "COA_biplot_by_region_all.pdf"),
       combined_geo, width = 18, height = 5, device = "pdf")

cat("Done.\n")
cat("Figures written to:", out_fig, "\n")
cat("Arrow correlations:", file.path(coa_dir, "COA_arrow_correlations.csv"), "\n")
cat("Axis contribution (wide):", file.path(coa_dir, "COA_axis_contribution_wide.csv"), "\n")
