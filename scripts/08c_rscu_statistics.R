# ============================================================
# 08c_rscu_statistics.R
# Purpose: Statistical tests for RSCU values.
#          For CDS, NSP and SP: Friedman test + Wilcoxon signed-rank
#          post-hoc + Benjamini-Hochberg (BH) correction.
#          For WEEV vs six hosts: one-sample Wilcoxon test + BH correction.
#          Note: CDS/NSP/SP are paired observations from the same sequence,
#                and each host provides a single pooled RSCU value, so
#                standard ANOVA is not applicable.
# Input:
#   results/tables/RSCU/WEEV/RSCU_CDS.csv
#   results/tables/RSCU/WEEV/RSCU_NSP.csv
#   results/tables/RSCU/WEEV/RSCU_SP.csv
#   results/tables/RSCU/RSCU.csv
# Output:
#   results/tables/RSCU/RSCU_stats_CDS_NSP_SP.csv
#   results/tables/RSCU/RSCU_stats_WEEV_hosts.csv
#   results/tables/RSCU/RSCU_stats_summary.csv
# Dependencies: here, dplyr, tidyr, readr, tibble, rstatix, purrr
# Random seed: 12345 (all tests are deterministic; seed kept per project convention)
# Run: Rscript scripts/08c_rscu_statistics.R
# ============================================================

set.seed(12345)  # random seed: 12345

suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr)
  library(readr); library(tibble); library(rstatix); library(purrr)
})

# ---------- 1. Read data ----------
read_wide <- function(path) {
  x <- read_csv(path, show_col_types = FALSE,
                col_types = cols(.default = col_character()))
  meta_cols <- intersect(c("Accession","AA","Codon"), names(x))
  num_cols  <- setdiff(names(x), meta_cols)
  x[num_cols] <- lapply(x[num_cols], as.numeric)
  x
}

wide_cds <- read_wide(here("results","tables","RSCU","WEEV","RSCU_CDS.csv"))
wide_nsp <- read_wide(here("results","tables","RSCU","WEEV","RSCU_NSP.csv"))
wide_sp  <- read_wide(here("results","tables","RSCU","WEEV","RSCU_SP.csv"))
rscu_all <- read_wide(here("results","tables","RSCU","RSCU.csv"))

# Exclude AUG and UGG
codons <- setdiff(names(wide_cds), "Accession")
codons <- codons[!codons %in% c("AUG","UGG")]

# Codon to amino acid
codon_to_aa <- c(
  UUU="Phe", UUC="Phe", UUA="Leu", UUG="Leu", UCU="Ser", UCC="Ser", UCA="Ser", UCG="Ser",
  UAU="Tyr", UAC="Tyr", UGU="Cys", UGC="Cys", UGG="Trp",
  CUU="Leu", CUC="Leu", CUA="Leu", CUG="Leu", CCU="Pro", CCC="Pro", CCA="Pro", CCG="Pro",
  CAU="His", CAC="His", CAA="Gln", CAG="Gln", CGU="Arg", CGC="Arg", CGA="Arg", CGG="Arg",
  AUU="Ile", AUC="Ile", AUA="Ile", AUG="Met", ACU="Thr", ACC="Thr", ACA="Thr", ACG="Thr",
  AAU="Asn", AAC="Asn", AAA="Lys", AAG="Lys", AGU="Ser", AGC="Ser", AGA="Arg", AGG="Arg",
  GUU="Val", GUC="Val", GUA="Val", GUG="Val", GCU="Ala", GCC="Ala", GCA="Ala", GCG="Ala",
  GAU="Asp", GAC="Asp", GAA="Glu", GAG="Glu", GGU="Gly", GGC="Gly", GGA="Gly", GGG="Gly"
)

# ---------- 2. CDS vs NSP vs SP ----------
res_i <- map_dfr(codons, function(cd) {
  df <- tibble(
    Accession = wide_cds$Accession,
    CDS = wide_cds[[cd]],
    NSP = wide_nsp[[cd]],
    SP  = wide_sp[[cd]]
  ) %>%
    pivot_longer(cols = c(CDS, NSP, SP), names_to = "Region", values_to = "RSCU") %>%
    mutate(Region = factor(Region, levels = c("CDS","NSP","SP")))

  fried <- friedman_test(df, RSCU ~ Region | Accession)

  pairs <- list(c("CDS","NSP"), c("CDS","SP"), c("NSP","SP"))
  post <- map_dfr(pairs, function(p) {
    sub <- df %>% filter(Region %in% p)
    w <- wilcox_test(sub, RSCU ~ Region, paired = TRUE)
    tibble(
      comparison = paste(p[1], "vs", p[2]),
      p = w$p
    )
  }) %>%
    pivot_wider(names_from = comparison, values_from = p,
                names_glue = "{comparison}_p")

  tibble(
    Codon = cd,
    AA = codon_to_aa[cd],
    Friedman_p = fried$p
  ) %>% bind_cols(post)
})

p_cols_i <- grep("_p$", names(res_i), value = TRUE)
res_i <- res_i %>%
  mutate(across(all_of(p_cols_i), ~ p.adjust(.x, method = "BH"),
                .names = "{.col}_adj"))

# ---------- 3. WEEV vs hosts ----------
host_cols <- c("Aedes aegypti","Culex tarsalis","Haemorhous mexicanus",
               "Passer domesticus","Equus caballus","Homo sapiens")

host_rscu <- rscu_all %>%
  filter(!Codon %in% c("AUG","UGG")) %>%
  select(Codon, all_of(host_cols))

res_ii <- map_dfr(host_cols, function(h) {
  map_dfr(codons, function(cd) {
    weev_vals <- wide_cds[[cd]]
    host_val  <- host_rscu[[h]][host_rscu$Codon == cd]
    if (length(host_val) != 1 || is.na(host_val)) return(NULL)
    w <- wilcox.test(weev_vals, mu = host_val)
    tibble(
      Codon = cd,
      AA = codon_to_aa[cd],
      Host = h,
      Wilcoxon_p = w$p.value,
      WEEV_median = median(weev_vals, na.rm = TRUE),
      Host_RSCU = host_val
    )
  })
})

res_ii <- res_ii %>%
  mutate(Wilcoxon_q = p.adjust(Wilcoxon_p, method = "BH"))

# ---------- 4. Summary ----------
sum_i <- res_i %>%
  summarise(
    comparison = "CDS vs NSP vs SP",
    n_tested = n(),
    n_sig_Friedman = sum(Friedman_p_adj < 0.05, na.rm = TRUE),
    sig_Friedman = paste(Codon[Friedman_p_adj < 0.05], collapse = ", ")
  )

sum_ii <- res_ii %>%
  group_by(Host) %>%
  summarise(
    n_tested = n(),
    n_sig = sum(Wilcoxon_q < 0.05, na.rm = TRUE),
    sig_codons = paste(Codon[Wilcoxon_q < 0.05], collapse = ", "),
    .groups = "drop"
  ) %>%
  mutate(Host = factor(Host, levels = host_cols)) %>%
  arrange(Host) %>%
  mutate(Host = as.character(Host)) %>%
  rename(comparison = Host)

summary_tbl <- bind_rows(
  sum_i %>% rename(n_sig = n_sig_Friedman, sig_codons = sig_Friedman),
  sum_ii
)

# ---------- 5. Write output ----------
out_dir <- here("results","tables","RSCU")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(res_i, here(out_dir, "RSCU_stats_CDS_NSP_SP.csv"))
write_csv(res_ii, here(out_dir, "RSCU_stats_WEEV_hosts.csv"))
write_csv(summary_tbl, here(out_dir, "RSCU_stats_summary.csv"))

# ---------- 6. Report ----------
cat("\n===== RSCU statistics completed =====\n")
cat("CDS vs NSP vs SP: ", nrow(res_i), " codons\n")
cat("  Significant by Friedman (q<0.05): ", sum(res_i$Friedman_p_adj < 0.05, na.rm = TRUE), "\n")
cat("WEEV vs hosts: ", nrow(res_ii), " rows (59 codons x 6 hosts)\n")
cat("  Significant (q<0.05): ", sum(res_ii$Wilcoxon_q < 0.05, na.rm = TRUE), "\n")
cat("\nOutput:\n")
cat("  results/tables/RSCU/RSCU_stats_CDS_NSP_SP.csv\n")
cat("  results/tables/RSCU/RSCU_stats_WEEV_hosts.csv\n")
cat("  results/tables/RSCU/RSCU_stats_summary.csv\n")