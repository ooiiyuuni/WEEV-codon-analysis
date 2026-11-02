# =============================================================================
# 09_sid_calculate.R
#
# Purpose : Calculate the codon usage dissimilarity index (SiD) between each
#           WEEV sequence (complete CDS, NSP, SP) and each of the six host
#           reference codon usage profiles.
#
#           SiD is defined as D(A,B) = (1 - R(A,B)) / 2, where R(A,B) is the
#           cosine similarity between the 59-dimensional RSCU vectors of the
#           viral sequence (A) and the host (B). Because all RSCU values are
#           non-negative, R is bounded in [0, 1] and D in [0, 0.5]. Higher D
#           indicates greater dissimilarity from the host reference.
#
# Input   : results/tables/RSCU/WEEV/RSCU_long_CDS.csv
#           results/tables/RSCU/WEEV/RSCU_long_NSP.csv
#           results/tables/RSCU/WEEV/RSCU_long_SP.csv
#           results/tables/RSCU/hosts/Aedes_aegypti_codon_freq.txt
#           results/tables/RSCU/hosts/Culex_tarsalis_codon_freq.txt
#           results/tables/RSCU/hosts/Haemorhous_mexicanus_codon_freq.txt
#           results/tables/RSCU/hosts/Passer_domesticus_codon_freq.txt
#           results/tables/RSCU/hosts/Equus_caballus_codon_freq.txt
#           results/tables/RSCU/hosts/Homo_sapiens_codon_freq.txt
#
# Output  : results/tables/SiD/SiD_long.csv
#           Columns: Accession, Region, Host, R, SiD
#           52 sequences x 3 regions x 6 hosts = 936 rows
#
# Dependencies: here, tidyverse
#
# Random seed : set.seed(12345) (no stochastic step; set for reproducibility
#               convention across all scripts in this project)
#
# Run     : Rscript scripts/09_sid_calculate.R
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
})

set.seed(12345)

# ---------------------------------------------------------------------------
# 1. Read WEEV per-sequence RSCU (long format) for three regions
# ---------------------------------------------------------------------------

read_weev_long <- function(path, region) {
  read_csv(path, show_col_types = FALSE) %>%
    select(Accession, Codon, RSCU) %>%
    filter(!Codon %in% c("AUG", "UGG")) %>%
    mutate(Region = region)
}

weev <- bind_rows(
  read_weev_long(here("results/tables/RSCU/WEEV/RSCU_long_CDS.csv"), "CDS"),
  read_weev_long(here("results/tables/RSCU/WEEV/RSCU_long_NSP.csv"), "NSP"),
  read_weev_long(here("results/tables/RSCU/WEEV/RSCU_long_SP.csv"),  "SP")
)

stopifnot(n_distinct(weev$Codon) == 59)
stopifnot(n_distinct(weev$Accession) == 52)

# ---------------------------------------------------------------------------
# 2. Read host codon frequencies and convert to RSCU
#    RSCU_ij = x_ij * n_i / sum_j(x_ij)
#    where n_i is the number of synonymous codons for amino acid i.
#    Using frequencies or counts gives identical RSCU.
# ---------------------------------------------------------------------------

host_files <- c(
  "Aedes aegypti"        = "Aedes_aegypti_codon_freq.txt",
  "Culex tarsalis"       = "Culex_tarsalis_codon_freq.txt",
  "Haemorhous mexicanus" = "Haemorhous_mexicanus_codon_freq.txt",
  "Passer domesticus"    = "Passer_domesticus_codon_freq.txt",
  "Equus caballus"       = "Equus_caballus_codon_freq.txt",
  "Homo sapiens"         = "Homo_sapiens_codon_freq.txt"
)

# Amino acid assignment for each sense codon (standard genetic code)
codon_aa <- c(
  UUU="Phe", UUC="Phe", UUA="Leu", UUG="Leu",
  CUU="Leu", CUC="Leu", CUA="Leu", CUG="Leu",
  AUU="Ile", AUC="Ile", AUA="Ile", AUG="Met",
  GUU="Val", GUC="Val", GUA="Val", GUG="Val",
  UCU="Ser", UCC="Ser", UCA="Ser", UCG="Ser",
  CCU="Pro", CCC="Pro", CCA="Pro", CCG="Pro",
  ACU="Thr", ACC="Thr", ACA="Thr", ACG="Thr",
  GCU="Ala", GCC="Ala", GCA="Ala", GCG="Ala",
  UAU="Tyr", UAC="Tyr", UAA="Stop", UAG="Stop",
  CAU="His", CAC="His", CAA="Gln", CAG="Gln",
  AAU="Asn", AAC="Asn", AAA="Lys", AAG="Lys",
  GAU="Asp", GAC="Asp", GAA="Glu", GAG="Glu",
  UGU="Cys", UGC="Cys", UGA="Stop", UGG="Trp",
  CGU="Arg", CGC="Arg", CGA="Arg", CGG="Arg",
  AGU="Ser", AGC="Ser", AGA="Arg", AGG="Arg",
  GGU="Gly", GGC="Gly", GGA="Gly", GGG="Gly"
)

read_host_rscu <- function(fname, host_name) {
  df <- read_tsv(
    here("results/tables/RSCU/hosts", fname),
    col_names = c("Codon", "Freq"),
    show_col_types = FALSE
  ) %>%
    filter(!Codon %in% c("AUG", "UGG")) %>%
    mutate(AA = codon_aa[Codon]) %>%
    group_by(AA) %>%
    mutate(n_syn    = n(),
           RSCU_host = Freq * n_syn / sum(Freq)) %>%
    ungroup() %>%
    select(Codon, RSCU_host) %>%
    mutate(Host = host_name)
  df
}

host_rscu <- purrr::map2_dfr(
  host_files, names(host_files), read_host_rscu
)

stopifnot(n_distinct(host_rscu$Codon) == 59)
stopifnot(n_distinct(host_rscu$Host) == 6)

# ---------------------------------------------------------------------------
# 3. Align by Codon and compute cosine similarity R and SiD = (1 - R) / 2
# ---------------------------------------------------------------------------

cosine_sim <- function(a, b) {
  sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
}

weev_wide <- weev %>%
  select(Accession, Region, Codon, RSCU) %>%
  tidyr::pivot_wider(names_from = Codon, values_from = RSCU)

host_wide <- host_rscu %>%
  select(Host, Codon, RSCU_host) %>%
  tidyr::pivot_wider(names_from = Codon, values_from = RSCU_host)

codon_order <- sort(unique(weev$Codon))
stopifnot(setequal(codon_order, sort(unique(host_rscu$Codon))))

weev_mat <- as.matrix(weev_wide[, codon_order])
host_mat <- as.matrix(host_wide[, codon_order])
rownames(weev_mat) <- paste(weev_wide$Accession, weev_wide$Region, sep = "__")
rownames(host_mat) <- host_wide$Host

R_mat <- outer(
  seq_len(nrow(weev_mat)),
  seq_len(nrow(host_mat)),
  Vectorize(function(i, j) cosine_sim(weev_mat[i, ], host_mat[j, ]))
)
rownames(R_mat) <- rownames(weev_mat)
colnames(R_mat) <- rownames(host_mat)

R_long <- as.data.frame(as.table(R_mat)) %>%
  rename(Vector = Var1, Host = Var2, R = Freq) %>%
  tidyr::separate(Vector, into = c("Accession", "Region"), sep = "__") %>%
  mutate(SiD = (1 - R) / 2) %>%
  select(Accession, Region, Host, R, SiD)

# ---------------------------------------------------------------------------
# 4. Sanity checks and output
# ---------------------------------------------------------------------------

stopifnot(nrow(R_long) == 52 * 3 * 6)
stopifnot(all(R_long$R   >= 0 & R_long$R   <= 1))
stopifnot(all(R_long$SiD >= 0 & R_long$SiD <= 0.5))

dir.create(here("results/tables/SiD"), recursive = TRUE, showWarnings = FALSE)
write_csv(R_long, here("results/tables/SiD/SiD_long.csv"))

message("SiD_long.csv written: ", nrow(R_long), " rows")
message("R range  : ", round(min(R_long$R),   4), " - ", round(max(R_long$R),   4))
message("SiD range: ", round(min(R_long$SiD), 4), " - ", round(max(R_long$SiD), 4))