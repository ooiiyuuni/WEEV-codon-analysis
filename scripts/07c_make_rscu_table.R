# ============================================================
# 07c_make_rscu_table.R
# Purpose: Merge WEEV RSCU and six host RSCU values into a single
#          AA/Codon/RSCU table. RSCU values are rounded to 4 decimals
#          (1 is written as 1.0000). Codons are ordered by amino acid group,
#          following the standard codon order used in the manuscript.
# Input:
#   results/tables/RSCU/WEEV/RSCU_long_CDS.csv
#   results/tables/RSCU/WEEV/RSCU_long_NSP.csv
#   results/tables/RSCU/WEEV/RSCU_long_SP.csv
#   results/tables/RSCU/hosts/*_codon_freq.txt
# Output:
#   results/tables/RSCU/RSCU.csv
# Dependencies: here, dplyr, tidyr, readr, tibble
# Random seed: 12345 (no stochastic step; kept per project convention)
# Run: Rscript scripts/07c_make_rscu_table.R
# ============================================================

set.seed(12345)  # random seed: 12345

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
})

# ---------- 1. Codon order (as in manuscript) ----------
codon_order <- c(
  "UUU","UUC",                                  # Phe
  "UUA","UUG","CUU","CUC","CUA","CUG",          # Leu
  "UCU","UCC","UCA","UCG","AGU","AGC",          # Ser
  "UAU","UAC",                                  # Tyr
  "UGU","UGC",                                  # Cys
  "CCU","CCC","CCA","CCG",                      # Pro
  "CAU","CAC",                                  # His
  "CAA","CAG",                                  # Gln
  "CGU","CGC","CGA","CGG","AGA","AGG",          # Arg
  "AUU","AUC","AUA",                            # Ile
  "AUG",                                        # Met
  "ACU","ACC","ACA","ACG",                      # Thr
  "AAU","AAC",                                  # Asn
  "AAA","AAG",                                  # Lys
  "GUU","GUC","GUA","GUG",                      # Val
  "GCU","GCC","GCA","GCG",                      # Ala
  "GAU","GAC",                                  # Asp
  "GAA","GAG",                                  # Glu
  "GGU","GGC","GGA","GGG",                      # Gly
  "UGG"                                         # Trp
)

codon_to_aa <- c(
  UUU="Phe", UUC="Phe", UUA="Leu", UUG="Leu",
  UCU="Ser", UCC="Ser", UCA="Ser", UCG="Ser",
  UAU="Tyr", UAC="Tyr", UGU="Cys", UGC="Cys", UGG="Trp",
  CUU="Leu", CUC="Leu", CUA="Leu", CUG="Leu",
  CCU="Pro", CCC="Pro", CCA="Pro", CCG="Pro",
  CAU="His", CAC="His", CAA="Gln", CAG="Gln",
  CGU="Arg", CGC="Arg", CGA="Arg", CGG="Arg",
  AUU="Ile", AUC="Ile", AUA="Ile", AUG="Met",
  ACU="Thr", ACC="Thr", ACA="Thr", ACG="Thr",
  AAU="Asn", AAC="Asn", AAA="Lys", AAG="Lys",
  AGU="Ser", AGC="Ser", AGA="Arg", AGG="Arg",
  GUU="Val", GUC="Val", GUA="Val", GUG="Val",
  GCU="Ala", GCC="Ala", GCA="Ala", GCG="Ala",
  GAU="Asp", GAC="Asp", GAA="Glu", GAG="Glu",
  GGU="Gly", GGC="Gly", GGA="Gly", GGG="Gly"
)

# ---------- 2. WEEV mean RSCU per codon ----------
weev_mean <- function(long_path, label) {
  read_csv(long_path, show_col_types = FALSE) %>%
    group_by(Codon) %>%
    summarise(!!label := mean(RSCU), .groups = "drop")
}

weev_cds <- weev_mean(here("results","tables","RSCU","WEEV","RSCU_long_CDS.csv"), "Complete")
weev_nsp <- weev_mean(here("results","tables","RSCU","WEEV","RSCU_long_NSP.csv"), "NSP")
weev_sp  <- weev_mean(here("results","tables","RSCU","WEEV","RSCU_long_SP.csv"),  "SP")

# ---------- 3. Host frequency to RSCU ----------
freq_to_rscu <- function(freq_path, label) {
  df <- read_tsv(freq_path, col_names = c("Codon","Freq"), show_col_types = FALSE)
  df$AA <- codon_to_aa[df$Codon]
  df %>%
    group_by(AA) %>%
    mutate(ngrp = n(),
           RSCU = if (first(ngrp) == 1) 1 else Freq / (sum(Freq) / first(ngrp))) %>%
    ungroup() %>%
    select(Codon, !!label := RSCU)
}

host_files <- tribble(
  ~label,                  ~file,
  "Aedes aegypti",         "Aedes_aegypti_codon_freq.txt",
  "Culex tarsalis",        "Culex_tarsalis_codon_freq.txt",
  "Haemorhous mexicanus",  "Haemorhous_mexicanus_codon_freq.txt",
  "Passer domesticus",     "Passer_domesticus_codon_freq.txt",
  "Equus caballus",        "Equus_caballus_codon_freq.txt",
  "Homo sapiens",          "Homo_sapiens_codon_freq.txt"
)

host_list <- lapply(seq_len(nrow(host_files)), function(i) {
  freq_to_rscu(
    here("results","tables","RSCU","hosts", host_files$file[i]),
    host_files$label[i]
  )
})

# ---------- 4. Merge ----------
base_tbl <- tibble(Codon = codon_order, AA = codon_to_aa[codon_order])

out_tbl <- base_tbl %>%
  left_join(weev_cds, by = "Codon") %>%
  left_join(weev_nsp, by = "Codon") %>%
  left_join(weev_sp,  by = "Codon")

for (h in host_list) out_tbl <- left_join(out_tbl, h, by = "Codon")

out_tbl <- out_tbl %>%
  select(AA, Codon, Complete, NSP, SP,
         `Aedes aegypti`, `Culex tarsalis`, `Haemorhous mexicanus`,
         `Passer domesticus`, `Equus caballus`, `Homo sapiens`)

# ---------- 5. Write output (4 decimals, 1 as 1.0000) ----------
out_dir <- here("results","tables","RSCU")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_fmt <- out_tbl %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))

write_csv(out_fmt, here(out_dir, "RSCU.csv"))

# ---------- 6. Report ----------
cat("\n===== RSCU table completed =====\n")
cat("Rows (expected 61): ", nrow(out_fmt), "\n")
cat("Columns (expected 11): ", ncol(out_fmt), "\n")
cat("\nOutput: results/tables/RSCU/RSCU.csv\n")