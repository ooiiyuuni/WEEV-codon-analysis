# =============================================================================
# 09c_sid_statistics.R
#
# Purpose : Test whether SiD differs (i) among the six host references and
#           (ii) among the three genomic regions (CDS, NSP, SP).
#           Host comparison: each WEEV sequence is compared against all six
#           hosts, so observations are paired across hosts; Friedman omnibus
#           followed by pairwise Wilcoxon signed-rank tests with
#           Benjamini-Hochberg correction.
#           Region comparison: for each host separately, the three regions are
#           paired within each sequence; Friedman omnibus followed by pairwise
#           Wilcoxon signed-rank tests with Benjamini-Hochberg correction.
#
# Input   : results/tables/SiD/SiD_long.csv
#           (columns: Accession, Region, Host, R, SiD)
#
# Output  : results/tables/SiD/SiD_stats.csv
#           Columns: Analysis, Group, Region/Host_scope, Test, Comparison,
#                    Statistic, df, P_value, BH_q
#
# Dependencies: here, tidyverse
#
# Random seed : set.seed(12345) (no stochastic step; project convention)
#
# Run     : Rscript scripts/09c_sid_statistics.R
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
})

set.seed(12345)

# ---------------------------------------------------------------------------
# 1. Read data and fix orders
# ---------------------------------------------------------------------------

sid <- read_csv(here("results/tables/SiD/SiD_long.csv"), show_col_types = FALSE)

host_order <- c(
  "Aedes aegypti",
  "Culex tarsalis",
  "Haemorhous mexicanus",
  "Passer domesticus",
  "Equus caballus",
  "Homo sapiens"
)

region_order <- c("CDS", "NSP", "SP")

sid <- sid %>%
  mutate(
    Host   = factor(Host,   levels = host_order),
    Region = factor(Region, levels = region_order)
  )

stopifnot(nrow(sid) == 52 * 3 * 6)

# ===========================================================================
# PART A. Host comparison (within each region)
# ===========================================================================

# ---------------------------------------------------------------------------
# A1. Friedman omnibus across hosts, per region
# ---------------------------------------------------------------------------

friedman_host <- list()

for (rg in region_order) {
  dat <- sid %>%
    filter(Region == rg) %>%
    select(Accession, Host, SiD) %>%
    pivot_wider(names_from = Host, values_from = SiD) %>%
    arrange(Accession)

  mat <- as.matrix(dat[, host_order])
  rownames(mat) <- dat$Accession

  ft <- friedman.test(mat)

  friedman_host[[rg]] <- tibble(
    Analysis   = "Host comparison",
    Group      = rg,
    Scope      = "All six hosts",
    Test       = "Friedman",
    Comparison = "All six hosts",
    Statistic  = unname(ft$statistic),
    df         = unname(ft$parameter),
    P_value    = ft$p.value,
    BH_q       = NA_real_
  )
}

friedman_host_df <- bind_rows(friedman_host)

# ---------------------------------------------------------------------------
# A2. Pairwise Wilcoxon across hosts, per region, BH within region
# ---------------------------------------------------------------------------

pairwise_host <- list()

for (rg in region_order) {
  wide <- sid %>%
    filter(Region == rg) %>%
    select(Accession, Host, SiD) %>%
    pivot_wider(names_from = Host, values_from = SiD) %>%
    arrange(Accession)

  pairs <- combn(host_order, 2, simplify = FALSE)

  rows <- lapply(pairs, function(pr) {
    wt <- wilcox.test(wide[[pr[1]]], wide[[pr[2]]],
                      paired = TRUE, exact = FALSE)
    tibble(
      Analysis   = "Host comparison",
      Group      = rg,
      Scope      = "All six hosts",
      Test       = "Wilcoxon signed-rank (paired)",
      Comparison = paste0(pr[1], " vs ", pr[2]),
      Statistic  = unname(wt$statistic),
      df         = NA_real_,
      P_value    = wt$p.value
    )
  })

  pairwise_host[[rg]] <- bind_rows(rows) %>%
    mutate(BH_q = p.adjust(P_value, method = "BH"))
}

pairwise_host_df <- bind_rows(pairwise_host)

# ===========================================================================
# PART B. Region comparison (within each host)
# ===========================================================================

# ---------------------------------------------------------------------------
# B1. Friedman omnibus across regions, per host
# ---------------------------------------------------------------------------

friedman_region <- list()

for (hs in host_order) {
  dat <- sid %>%
    filter(Host == hs) %>%
    select(Accession, Region, SiD) %>%
    pivot_wider(names_from = Region, values_from = SiD) %>%
    arrange(Accession)

  mat <- as.matrix(dat[, region_order])
  rownames(mat) <- dat$Accession

  ft <- friedman.test(mat)

  friedman_region[[hs]] <- tibble(
    Analysis   = "Region comparison",
    Group      = hs,
    Scope      = "All three regions",
    Test       = "Friedman",
    Comparison = "All three regions",
    Statistic  = unname(ft$statistic),
    df         = unname(ft$parameter),
    P_value    = ft$p.value,
    BH_q       = NA_real_
  )
}

friedman_region_df <- bind_rows(friedman_region)

# ---------------------------------------------------------------------------
# B2. Pairwise Wilcoxon across regions, per host, BH within host
# ---------------------------------------------------------------------------

pairwise_region <- list()

for (hs in host_order) {
  wide <- sid %>%
    filter(Host == hs) %>%
    select(Accession, Region, SiD) %>%
    pivot_wider(names_from = Region, values_from = SiD) %>%
    arrange(Accession)

  pairs <- combn(region_order, 2, simplify = FALSE)

  rows <- lapply(pairs, function(pr) {
    wt <- wilcox.test(wide[[pr[1]]], wide[[pr[2]]],
                      paired = TRUE, exact = FALSE)
    tibble(
      Analysis   = "Region comparison",
      Group      = hs,
      Scope      = "All three regions",
      Test       = "Wilcoxon signed-rank (paired)",
      Comparison = paste0(pr[1], " vs ", pr[2]),
      Statistic  = unname(wt$statistic),
      df         = NA_real_,
      P_value    = wt$p.value
    )
  })

  pairwise_region[[hs]] <- bind_rows(rows) %>%
    mutate(BH_q = p.adjust(P_value, method = "BH"))
}

pairwise_region_df <- bind_rows(pairwise_region)

# ===========================================================================
# PART C. Combine and write
# ===========================================================================

stats_out <- bind_rows(
  friedman_host_df,
  pairwise_host_df,
  friedman_region_df,
  pairwise_region_df
) %>%
  mutate(
    Statistic = round(Statistic, 4),
    df        = ifelse(is.na(df), NA, df),
    P_value   = signif(P_value, 4),
    BH_q      = ifelse(is.na(BH_q), NA, signif(BH_q, 4))
  ) %>%
  select(Analysis, Group, Scope, Test, Comparison, Statistic, df, P_value, BH_q)

dir.create(here("results/tables/SiD"), recursive = TRUE, showWarnings = FALSE)
write_csv(stats_out, here("results/tables/SiD/SiD_stats.csv"))

message("SiD_stats.csv written: ", nrow(stats_out), " rows")

# ===========================================================================
# PART D. Console summaries
# ===========================================================================

cat("\n=== Friedman: hosts within each region ===\n")
print(friedman_host_df %>% select(Group, Statistic, df, P_value))

cat("\n=== Friedman: regions within each host ===\n")
print(friedman_region_df %>% select(Group, Statistic, df, P_value))

cat("\n=== Pairwise (hosts) with BH q >= 0.05 (if any) ===\n")
print(pairwise_host_df %>% filter(BH_q >= 0.05) %>%
        select(Group, Comparison, Statistic, P_value, BH_q))

cat("\n=== Pairwise (regions) with BH q >= 0.05 (if any) ===\n")
print(pairwise_region_df %>% filter(BH_q >= 0.05) %>%
        select(Group, Comparison, Statistic, P_value, BH_q))