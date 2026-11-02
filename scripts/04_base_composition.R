# ============================================================
# 04_base_composition.R
# Function: Calculate base composition (A,U,G,C,GC,AUGC,Deviation)
# Input: data/WEEV/{cds,nsp,sp}.fasta
# Output: results/tables/BaseComposition_{CDS,NSP,SP}.csv
# Dependencies: here, seqinr, tidyverse
# Run: Rscript scripts/04_base_composition.R
# ============================================================

library(here)
library(seqinr)
library(tidyverse)

read_weev <- function(fname, suffix) {
  s <- read.fasta(here("data", "WEEV", fname), seqtype = "DNA", as.string = TRUE)
  names(s) <- sub(suffix, "", names(s))
  s
}

cds <- read_weev("cds.fasta", "_complete_CDS$")
nsp <- read_weev("nsp.fasta", "_nsp$")
sp  <- read_weev("sp.fasta",  "_sp$")

base_composition <- function(seq_list) {
  tibble(
    Accession = names(seq_list),
    Length    = sapply(seq_list, function(x) nchar(x[1])),
    A = sapply(seq_list, function(x) str_count(toupper(x[1]), "A")),
    U = sapply(seq_list, function(x) str_count(toupper(x[1]), "T")),
    G = sapply(seq_list, function(x) str_count(toupper(x[1]), "G")),
    C = sapply(seq_list, function(x) str_count(toupper(x[1]), "C"))
  ) %>%
    mutate(
      `A%`        = A / Length * 100,
      `U%`        = U / Length * 100,
      `G%`        = G / Length * 100,
      `C%`        = C / Length * 100,
      `GC%`       = (G + C) / Length * 100,
      `AUGC%`     = (A + U + G + C) / Length * 100,
      `Deviation` = 100 - `AUGC%`
    ) %>%
    select(Accession, Length, `A%`, `U%`, `G%`, `C%`, `GC%`, `AUGC%`, `Deviation`)
}

append_summary <- function(df, file_name, dataset_name) {
  mean_row <- df %>% summarise(across(where(is.numeric), ~ mean(.x))) %>%
    mutate(Accession = paste0(dataset_name, "_Mean"), .before = 1)
  sd_row <- df %>% summarise(across(where(is.numeric), ~ sd(.x))) %>%
    mutate(Accession = paste0(dataset_name, "_SD"), .before = 1)

  out <- bind_rows(df, mean_row, sd_row) %>%
    mutate(
      Length = ifelse(grepl("_Mean|_SD", Accession), sprintf("%.1f", Length), sprintf("%.0f", Length)),
      across(c(`A%`, `U%`, `G%`, `C%`, `GC%`, `AUGC%`, `Deviation`), ~ sprintf("%.4f", .x))
    )
  
  out_dir <- here("results", "tables")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(out, file.path(out_dir, file_name), row.names = FALSE, quote = FALSE)
}

append_summary(base_composition(cds), "BaseComposition_CDS.csv", "CDS")
append_summary(base_composition(nsp), "BaseComposition_NSP.csv", "NSP")
append_summary(base_composition(sp),  "BaseComposition_SP.csv",  "SP")

cat("Base composition saved to results/tables/\n")