# ============================================================
# 16_coa.R
# Purpose : Correspondence analysis (COA) of codon usage variation
#           across 52 WEEV sequences for complete CDS, NSP and SP,
#           based on RSCU values of 59 sense codons
#           (AUG and UGG excluded, as in the Methods).
# Input   : results/tables/RSCU/WEEV/RSCU_long_{CDS,NSP,SP}.csv
# Output  : results/tables/COA/COA_coords_{CDS,NSP,SP}.csv
#           results/tables/COA/COA_axis_contribution.csv
#           results/tables/COA/COA_codon_contribution_{CDS,NSP,SP}.csv
# Dependencies: here, tidyverse
# Random seed : not applicable (no stochastic step)
# Run     : Rscript scripts/16_coa.R
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
})

long_files <- c(
  CDS = here("results", "tables", "RSCU", "WEEV", "RSCU_long_CDS.csv"),
  NSP = here("results", "tables", "RSCU", "WEEV", "RSCU_long_NSP.csv"),
  SP  = here("results", "tables", "RSCU", "WEEV", "RSCU_long_SP.csv")
)

out_dir <- here("results", "tables", "COA")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

for (f in long_files) stopifnot(file.exists(f))

excluded_codons <- c("AUG", "UGG")

read_rscu_matrix <- function(path) {
  df <- read_csv(path, show_col_types = FALSE)
  stopifnot(all(c("Accession", "Codon", "AA", "RSCU") %in% names(df)))
  stopifnot(!any(is.na(df$RSCU)))
  df <- df %>% filter(!Codon %in% excluded_codons)
  mat <- df %>%
    select(Accession, Codon, RSCU) %>%
    pivot_wider(names_from = Codon, values_from = RSCU) %>%
    column_to_rownames("Accession") %>%
    as.matrix()
  mat
}

mats <- lapply(long_files, read_rscu_matrix)

for (nm in names(mats)) {
  stopifnot(nrow(mats[[nm]]) == 52)
  stopifnot(ncol(mats[[nm]]) == 59)
  stopifnot(!any(is.na(mats[[nm]])))
}

stopifnot(identical(rownames(mats$CDS), rownames(mats$NSP)))
stopifnot(identical(rownames(mats$CDS), rownames(mats$SP)))
stopifnot(identical(colnames(mats$CDS), colnames(mats$NSP)))
stopifnot(identical(colnames(mats$CDS), colnames(mats$SP)))

run_coa <- function(mat) {
  P <- mat / sum(mat)
  r <- rowSums(P)
  c <- colSums(P)
  S <- (P - outer(r, c)) / sqrt(outer(r, c))
  sv <- svd(S)
  eig <- sv$d^2
  contrib <- eig / sum(eig) * 100
  row_coords <- sv$u %*% diag(sv$d)
  rownames(row_coords) <- rownames(mat)
  colnames(row_coords) <- paste0("Axis", seq_len(ncol(row_coords)))
  col_coords <- sv$v %*% diag(sv$d)
  rownames(col_coords) <- colnames(mat)
  colnames(col_coords) <- paste0("Axis", seq_len(ncol(col_coords)))
  list(
    row_coords = row_coords,
    col_coords = col_coords,
    contribution = contrib
  )
}

coa_results <- lapply(mats, run_coa)

for (region in c("CDS", "NSP", "SP")) {
  rc <- coa_results[[region]]$row_coords
  out <- tibble(Accession = rownames(rc), Axis1 = rc[, 1], Axis2 = rc[, 2])
  write_csv(out, file.path(out_dir, paste0("COA_coords_", region, ".csv")))
}

max_axis <- 10
contrib_tab <- map_dfr(c("CDS", "NSP", "SP"), function(region) {
  v <- coa_results[[region]]$contribution
  k <- min(max_axis, length(v))
  tibble(
    Region = region,
    Axis = seq_len(k),
    Contribution_percent = v[seq_len(k)],
    Cumulative_percent = cumsum(v[seq_len(k)])
  )
})
write_csv(contrib_tab, file.path(out_dir, "COA_axis_contribution.csv"))

for (region in c("CDS", "NSP", "SP")) {
  cc <- coa_results[[region]]$col_coords
  out <- tibble(
    Codon = rownames(cc),
    Axis1 = cc[, 1],
    Axis2 = cc[, 2]
  )
  write_csv(out, file.path(out_dir, paste0("COA_codon_contribution_", region, ".csv")))
}

cat("Done.\n")
cat("Files written to:", out_dir, "\n")
