\# ============================================================
# 07b_host_codon_freq.R
# Purpose: Calculate codon usage frequency (61 sense codons) for six hosts
#          from locally stored CDS FASTA files. The denominator is the total
#          number of sense codons, so the 61 frequencies sum to 1.00.
#          Host CDS files are not included in the repository due to size.
#          Pre-computed outputs are provided in results/tables/RSCU/hosts/.
#          To rerun this script, download the original CDS files from NCBI
#          RefSeq (accessions listed in README) and update the paths below.
# Input (local paths, not uploaded):
#   ~/Desktop/WEEV/宿主/ajyw_cds.fna
#   ~/Desktop/Culex tarsalis提取cds/CtarK1_cds_v2.fa
#   ~/Desktop/WEEV/宿主/jzq_cds.fna
#   ~/Desktop/WEEV/宿主/jmq_cds.fna
#   ~/Desktop/WEEV/宿主/m_cds.fna
#   ~/Desktop/WEEV/宿主/homo_cds.fna
# Output:
#   results/tables/RSCU/hosts/Aedes_aegypti_codon_freq.txt
#   results/tables/RSCU/hosts/Culex_tarsalis_codon_freq.txt
#   results/tables/RSCU/hosts/Haemorhous_mexicanus_codon_freq.txt
#   results/tables/RSCU/hosts/Passer_domesticus_codon_freq.txt
#   results/tables/RSCU/hosts/Equus_caballus_codon_freq.txt
#   results/tables/RSCU/hosts/Homo_sapiens_codon_freq.txt
#   results/tables/RSCU/hosts/host_codon_freq_check.csv
# Dependencies: here, Biostrings, dplyr, readr, tibble
# Random seed: none (no stochastic step)
# Run: Rscript scripts/07b_host_codon_freq.R
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(Biostrings)
  library(dplyr)
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

# ---------- 2. Codon frequency for one FASTA ----------
calc_codon_freq <- function(fasta_path) {
  dna <- readDNAStringSet(fasta_path)
  total <- 0L
  counts <- setNames(integer(length(sense_codons)), sense_codons)
  for (i in seq_along(dna)) {
    s <- gsub("T", "U", as.character(dna[[i]]))
    n <- nchar(s)
    if (n < 3) next
    starts <- seq(1, n - 2, by = 3)
    codons <- substring(s, starts, starts + 2)
    codons <- codons[!codons %in% stop_codons]
    if (length(codons) == 0) next
    tab <- table(codons)
    hit <- intersect(names(tab), sense_codons)
    counts[hit] <- counts[hit] + as.integer(tab[hit])
    total <- total + length(codons)
  }
  tibble(Codon = sense_codons, Freq = as.numeric(counts / total))
}

# ---------- 3. Host CDS paths (local, not uploaded) ----------
hosts <- tribble(
  ~species,                 ~path,
  "Aedes_aegypti",          "~/Desktop/WEEV/宿主/ajyw_cds.fna",
  "Culex_tarsalis",         "~/Desktop/Culex tarsalis提取cds/CtarK1_cds_v2.fa",
  "Haemorhous_mexicanus",   "~/Desktop/WEEV/宿主/jzq_cds.fna",
  "Passer_domesticus",      "~/Desktop/WEEV/宿主/jmq_cds.fna",
  "Equus_caballus",         "~/Desktop/WEEV/宿主/m_cds.fna",
  "Homo_sapiens",           "~/Desktop/WEEV/宿主/homo_cds.fna"
)

# ---------- 4. Main loop ----------
out_dir <- here("results","tables","RSCU","hosts")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

qc <- list()
for (i in seq_len(nrow(hosts))) {
  sp   <- hosts$species[i]
  path <- path.expand(hosts$path[i])
  if (!file.exists(path)) stop("File not found: ", path)
  message("Processing: ", sp, " (", basename(path), ")")
  freq <- calc_codon_freq(path)
  write_lines(sprintf("%s\t%.8f", freq$Codon, freq$Freq),
              here(out_dir, paste0(sp, "_codon_freq.txt")))
  qc[[sp]] <- tibble(
    species    = sp,
    input_file = basename(path),
    n_rows     = nrow(freq),
    freq_sum   = sum(freq$Freq),
    min_freq   = min(freq$Freq),
    max_freq   = max(freq$Freq)
  )
}

qc_df <- bind_rows(qc)
write_csv(qc_df, here(out_dir, "host_codon_freq_check.csv"))

# ---------- 5. Report ----------
cat("\n===== Host codon frequency completed =====\n")
print(qc_df)
cat("\nOutput directory: results/tables/RSCU/hosts/\n")