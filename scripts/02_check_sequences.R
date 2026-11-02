# ============================================================
# 02_check_sequences.R
# Description: Validate that CDS, NSP, and SP sequences correspond one-to-one,
#              and check length relation: CDS = NSP + SP
# Input:  data/WEEV/cds.fasta, data/WEEV/nsp.fasta, data/WEEV/sp.fasta
# Output: results/tables/sequence_length_check.csv
# Dependencies: here, seqinr, tidyverse
# Random seed: Not applicable
# Run: Rscript scripts/02_check_sequences.R
# ============================================================

library(here)
library(seqinr)
library(tidyverse)

cds <- read.fasta(here("data", "WEEV", "cds.fasta"), seqtype = "DNA", as.string = TRUE)
nsp <- read.fasta(here("data", "WEEV", "nsp.fasta"), seqtype = "DNA", as.string = TRUE)
sp  <- read.fasta(here("data", "WEEV", "sp.fasta"),  seqtype = "DNA", as.string = TRUE)

names(cds) <- sub("_complete_CDS$", "", names(cds))
names(nsp) <- sub("_nsp$", "", names(nsp))
names(sp)  <- sub("_sp$", "", names(sp))

cat("CDS sequences:", length(cds), "\n")
cat("NSP sequences:", length(nsp), "\n")
cat("SP  sequences:", length(sp), "\n")

stopifnot(setequal(names(cds), names(nsp)))
stopifnot(setequal(names(cds), names(sp)))
cat("Accessions match across all three sets.\n")

ids <- sort(names(cds))

check <- tibble(
  accession          = ids,
  length_CDS         = sapply(ids, function(i) nchar(cds[[i]][1])),
  length_NSP         = sapply(ids, function(i) nchar(nsp[[i]][1])),
  length_SP          = sapply(ids, function(i) nchar(sp[[i]][1]))
) %>%
  mutate(
    length_NSP_plus_SP = length_NSP + length_SP,
    diff = length_CDS - length_NSP_plus_SP,
    ok   = diff == 0
  )

dir.create(here("results", "tables"), showWarnings = FALSE, recursive = TRUE)
write.csv(check, here("results", "tables", "sequence_length_check.csv"), row.names = FALSE)

cat("\n===== Sequence length check summary =====\n")
cat("CDS length range:", min(check$length_CDS), "-", max(check$length_CDS), "bp\n")
cat("NSP length range:", min(check$length_NSP), "-", max(check$length_NSP), "bp\n")
cat("SP  length range:", min(check$length_SP),  "-", max(check$length_SP),  "bp\n")
cat("CDS = NSP + SP for all sequences?", all(check$ok), "\n")
cat("Number of mismatches:", sum(!check$ok), "\n")
cat("\nResults saved to results/tables/sequence_length_check.csv\n")