# ============================================================
# 07_rscu_weev.R
# Purpose: Calculate RSCU for each WEEV CDS, NSP and SP sequence
#          from FASTA files. Output order follows the input FASTA.
# Input:
#   data/WEEV/cds.fasta
#   data/WEEV/nsp.fasta
#   data/WEEV/sp.fasta
# Output:
#   results/tables/RSCU/WEEV/RSCU_CDS.csv
#   results/tables/RSCU/WEEV/RSCU_NSP.csv
#   results/tables/RSCU/WEEV/RSCU_SP.csv
#   results/tables/RSCU/WEEV/RSCU_long_CDS.csv
#   results/tables/RSCU/WEEV/RSCU_long_NSP.csv
#   results/tables/RSCU/WEEV/RSCU_long_SP.csv
#   results/tables/RSCU/WEEV/RSCU_matrix_check.csv
# Dependencies: here, Biostrings, dplyr, tidyr, readr, tibble
# Random seed: 12345 (no stochastic step; kept per project convention)
# Run: Rscript scripts/07_rscu_weev.R
# ============================================================

set.seed(12345)  # random seed: 12345

suppressPackageStartupMessages({
  library(here)
  library(Biostrings)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
})

# ---------- 1. Codon definitions ----------
sense_codons <- c(
  "UUU","UUC","UUA","UUG","UCU","UCC","UCA","UCG",
  "UAU","UAC","UGU","UGC","UGG",
  "CUU","CUC","CUA","CUG","CCU","CCC","CCA","CCG",
  "CAU","CAC","CAA","CAG","CGU","CGC","CGA","CGG",
  "AUU","AUC","AUA","AUG","ACU","ACC","ACA","ACG",
  "AAU","AAC","AAA","AAG","AGU","AGC","AGA","AGG",
  "GUU","GUC","GUA","GUG","GCU","GCC","GCA","GCG",
  "GAU","GAC","GAA","GAG","GGU","GGC","GGA","GGG"
)
stop_codons <- c("UAA","UAG","UGA")

codon_to_aa <- c(
  UUU="Phe", UUC="Phe", UUA="Leu", UUG="Leu",
  UCU="Ser", UCC="Ser", UCA="Ser", UCG="Ser",
  UAU="Tyr", UAC="Tyr", UAA="Stop", UAG="Stop",
  UGU="Cys", UGC="Cys", UGA="Stop", UGG="Trp",
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

# ---------- 2. RSCU for a single sequence ----------
calc_rscu_one <- function(seq_str, id) {
  s <- gsub("T", "U", toupper(seq_str))
  n <- nchar(s)
  if (n < 3) return(NULL)
  starts <- seq(1, n - 2, by = 3)
  codons <- substring(s, starts, starts + 2)
  codons <- codons[!codons %in% stop_codons]
  if (length(codons) == 0) return(NULL)

  tab <- table(codons)
  counts <- setNames(integer(length(sense_codons)), sense_codons)
  hit <- intersect(names(tab), sense_codons)
  counts[hit] <- as.integer(tab[hit])

  aa_vec <- codon_to_aa[sense_codons]
  rscu <- numeric(length(sense_codons))
  names(rscu) <- sense_codons
  for (aa in unique(aa_vec)) {
    idx <- which(aa_vec == aa)
    if (length(idx) == 1) {
      rscu[idx] <- 1
    } else {
      grp_sum <- sum(counts[idx])
      rscu[idx] <- if (grp_sum == 0) 0 else counts[idx] / (grp_sum / length(idx))
    }
  }

  tibble(Accession = id, Codon = sense_codons, RSCU = as.numeric(rscu))
}

# ---------- 3. Batch processing (FASTA order preserved) ----------
calc_rscu_fasta <- function(fasta_path, suffix) {
  dna <- readDNAStringSet(fasta_path)
  ids <- sub(suffix, "", names(dna))
  res <- vector("list", length(dna))
  for (i in seq_along(dna)) {
    res[[i]] <- calc_rscu_one(as.character(dna[[i]]), ids[i])
  }
  bind_rows(res)
}

# ---------- 4. Main ----------
long_cds <- calc_rscu_fasta(here("data","WEEV","cds.fasta"), "_complete_CDS$")
long_nsp <- calc_rscu_fasta(here("data","WEEV","nsp.fasta"), "_nsp$")
long_sp  <- calc_rscu_fasta(here("data","WEEV","sp.fasta"),  "_sp$")

add_aa <- function(df) {
  df %>% mutate(AA = codon_to_aa[Codon]) %>% select(Accession, Codon, AA, RSCU)
}
long_cds <- add_aa(long_cds)
long_nsp <- add_aa(long_nsp)
long_sp  <- add_aa(long_sp)

to_wide <- function(long_df) {
  long_df %>%
    select(Accession, Codon, RSCU) %>%
    pivot_wider(names_from = Codon, values_from = RSCU)
}
wide_cds <- to_wide(long_cds)
wide_nsp <- to_wide(long_nsp)
wide_sp  <- to_wide(long_sp)

qc <- tibble(
  dataset  = c("CDS","NSP","SP"),
  n_rows   = c(nrow(wide_cds), nrow(wide_nsp), nrow(wide_sp)),
  n_codons = c(ncol(wide_cds) - 1, ncol(wide_nsp) - 1, ncol(wide_sp) - 1),
  n_NA     = c(sum(is.na(wide_cds)), sum(is.na(wide_nsp)), sum(is.na(wide_sp))),
  n_long   = c(nrow(long_cds), nrow(long_nsp), nrow(long_sp))
)

# ---------- 5. Write output ----------
out_dir <- here("results","tables","RSCU","WEEV")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(wide_cds, here(out_dir, "RSCU_CDS.csv"))
write_csv(wide_nsp, here(out_dir, "RSCU_NSP.csv"))
write_csv(wide_sp,  here(out_dir, "RSCU_SP.csv"))

write_csv(long_cds, here(out_dir, "RSCU_long_CDS.csv"))
write_csv(long_nsp, here(out_dir, "RSCU_long_NSP.csv"))
write_csv(long_sp,  here(out_dir, "RSCU_long_SP.csv"))

write_csv(qc, here(out_dir, "RSCU_matrix_check.csv"))

# ---------- 6. Report ----------
cat("\n===== WEEV RSCU matrix completed =====\n")
print(qc)
cat("\nOutput directory: results/tables/RSCU/WEEV/\n")