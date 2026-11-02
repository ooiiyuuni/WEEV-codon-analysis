、# ============================================================
# 03_degenerate_bases.R
# Function: Identify IUPAC degenerate bases in CDS sequences
# Input: data/WEEV/cds.fasta
# Output: results/tables/DegenerateBases.csv
# Dependencies: here, seqinr, tidyverse
# Run: Rscript scripts/03_degenerate_bases.R
# ============================================================

library(here)
library(seqinr)
library(tidyverse)

cds <- read.fasta(here("data", "WEEV", "cds.fasta"), seqtype = "DNA", as.string = TRUE)
names(cds) <- sub("_complete_CDS$", "", names(cds))

degenerate_table <- map_dfr(names(cds), function(id) {
  seq_str <- toupper(cds[[id]][1])
  chars   <- strsplit(seq_str, "")[[1]]
  pos     <- which(!chars %in% c("A", "C", "G", "T"))
  if (length(pos) == 0) return(NULL)
  bases <- chars[pos]
  
  tibble(
    Accession                 = id,
    `Length (bp)`             = nchar(seq_str),
    `No. of degenerate bases` = length(pos),
    `Degenerate bases`        = paste(sort(unique(bases)), collapse = "; "),
    Positions                 = paste(pos, collapse = "; "),
    `No. of affected codons`  = length(pos),
    `% of codons`             = round(length(pos) / (nchar(seq_str) / 3) * 100, 3)
  )
})

total_length_all <- sum(sapply(cds, function(x) nchar(x[1])))

total_row <- tibble(
  Accession                 = "Total",
  `Length (bp)`             = total_length_all,
  `No. of degenerate bases` = sum(degenerate_table$`No. of degenerate bases`),
  `Degenerate bases`        = "K; M; R; W; Y",
  Positions                 = "—",
  `No. of affected codons`  = sum(degenerate_table$`No. of affected codons`),
  `% of codons`             = round(sum(degenerate_table$`No. of affected codons`) / (total_length_all / 3) * 100, 3)
)

out <- bind_rows(degenerate_table, total_row)

dir.create(here("results", "tables"), showWarnings = FALSE, recursive = TRUE)
write.csv(out, here("results", "tables", "DegenerateBases.csv"), row.names = FALSE, quote = FALSE)

cat("Degenerate bases summary saved to results/tables/DegenerateBases.csv\n")