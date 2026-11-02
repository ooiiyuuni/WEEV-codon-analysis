# ============================================================
# Purpose : Prepare input files for CAIcal web server.
#           Remove codons with IUPAC bases and all stop codons.
# Input   : data/WEEV/{cds,nsp,sp}.fasta
#           results/tables/RSCU/hosts/*_codon_freq.txt
# Output  : results/tables/CAI/ (web_input fasta and host txt)
# Dependencies: here, seqinr, tidyverse
# Random seed  : no stochastic step
# Run     : Rscript scripts/13_prepare_web_inputs.R
# ============================================================

library(here)
library(seqinr)
library(tidyverse)

# ---- 0. Paths ----------------------------------------------
fasta_paths <- c(
  CDS = here("data", "WEEV", "cds.fasta"),
  NSP = here("data", "WEEV", "nsp.fasta"),
  SP  = here("data", "WEEV", "sp.fasta")
)
host_dir <- here("results", "tables", "RSCU", "hosts")
out_cai  <- here("results", "tables", "CAI")
dir.create(out_cai, showWarnings = FALSE, recursive = TRUE)

host_files <- c(
  "Aedes_aegypti"        = "Aedes_aegypti_codon_freq.txt",
  "Culex_tarsalis"       = "Culex_tarsalis_codon_freq.txt",
  "Haemorhous_mexicanus" = "Haemorhous_mexicanus_codon_freq.txt",
  "Passer_domesticus"    = "Passer_domesticus_codon_freq.txt",
  "Equus_caballus"       = "Equus_caballus_codon_freq.txt",
  "Homo_sapiens"         = "Homo_sapiens_codon_freq.txt"
)

# ---- 1. Clean WEEV fasta (直接删除) -------------------------
clean_fasta_delete <- function(in_file, out_file) {
  seqs <- read.fasta(in_file, as.string = TRUE, seqtype = "DNA", forceDNAtolower = FALSE)
  out <- list()
  for (i in seq_along(seqs)) {
    s <- toupper(as.character(seqs[[i]]))
    
    # 按 3 个碱基切分
    starts <- seq(1, nchar(s)-2, by=3)
    codons <- substring(s, starts, starts+2)
    
    # 1. 剔除含有非 ACGT 的密码子
    codons <- codons[grepl("^[ACGT]{3}$", codons)]
    
    # 2. 剔除所有终止密码子 (TAA, TAG, TGA)
    codons <- codons[!codons %in% c("TAA", "TAG", "TGA")]
    
    # 拼回去
    out[[i]] <- paste(codons, collapse="")
  }
  names(out) <- sub("_(complete_CDS|nsp|sp)$", "", names(seqs))
  write.fasta(out, names = names(out), file.out = out_file)
}

for (rg in names(fasta_paths)) {
  out_f <- file.path(out_cai, paste0("WEEV_", rg, "_clean.fasta"))
  clean_fasta_delete(fasta_paths[[rg]], out_f)
}

# ---- 2. Format host codon usage tables ---------------------
codon_order <- c("TTT","TTC","TTA","TTG","CTT","CTC","CTA","CTG","ATT","ATC","ATA","ATG","GTT","GTC","GTA","GTG","TCT","TCC","TCA","TCG","CCT","CCC","CCA","CCG","ACT","ACC","ACA","ACG","GCT","GCC","GCA","GCG","TAT","TAC","TAA","TAG","CAT","CAC","CAA","CAG","AAT","AAC","AAA","AAG","GAT","GAC","GAA","GAG","TGT","TGC","TGA","TGG","CGT","CGC","CGA","CGG","AGT","AGC","AGA","AGG","GGT","GGC","GGA","GGG")

for (h in names(host_files)) {
  d <- read.table(file.path(host_dir, host_files[h]), header=FALSE)
  d$V1 <- chartr("U","T",d$V1)
  
  # 补齐 3 个终止密码子，频率和计数全填 0
  stops <- data.frame(V1 = c("TAA", "TAG", "TGA"), V2 = c(0, 0, 0))
  d <- rbind(d, stops)
  d <- d[match(codon_order, d$V1), ]
  
  # 生成网页要求的格式
  total_codons <- 10000000
  lines <- paste0(d$V1, " ", sprintf("%.1f", d$V2 * 1000), "(", round(d$V2 * total_codons), ")")
  mat <- matrix(lines, ncol=4, byrow=TRUE)
  
  out_file <- file.path(out_cai, paste0("host_", h, ".txt"))
  write.table(mat, file=out_file, quote=FALSE, row.names=FALSE, col.names=FALSE)
}

message("Done. Web inputs saved to: ", out_cai)
