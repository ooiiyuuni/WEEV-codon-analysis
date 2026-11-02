# ============================================================
# 15_rad.R
# Purpose : Calculate relative abundance of dinucleotides (RAD)
#           for WEEV complete CDS, NSP and SP coding sequences,
#           and generate the summary table and barplot figure.
# Input   : data/WEEV/cds.fasta, data/WEEV/nsp.fasta, data/WEEV/sp.fasta
# Output  : results/tables/RAD/RAD_per_sequence_{CDS,NSP,SP}.csv
#           results/tables/RAD/RAD_summary.csv
#           results/figures/RAD_barplot.pdf
# Dependencies: here, tidyverse, Biostrings
# Random seed : not applicable (no stochastic step)
# Run     : Rscript scripts/15_rad.R
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(Biostrings)
})

fasta_files <- c(
  CDS = here("data", "WEEV", "cds.fasta"),
  NSP = here("data", "WEEV", "nsp.fasta"),
  SP  = here("data", "WEEV", "sp.fasta")
)

out_tab <- here("results", "tables", "RAD")
out_fig <- here("results", "figures")
if (!dir.exists(out_tab)) dir.create(out_tab, recursive = TRUE)
if (!dir.exists(out_fig)) dir.create(out_fig, recursive = TRUE)

for (f in fasta_files) {
  stopifnot(file.exists(f))
}

read_weev <- function(path) {
  dna <- readDNAStringSet(path)
  names(dna) <- sub("_(complete_CDS|nsp|sp)$", "", names(dna))
  dna
}

seqs <- lapply(fasta_files, read_weev)

for (nm in names(seqs)) {
  stopifnot(length(seqs[[nm]]) == 52)
  stopifnot(!any(is.na(names(seqs[[nm]]))))
  stopifnot(!any(duplicated(names(seqs[[nm]]))))
}

stopifnot(identical(names(seqs$CDS), names(seqs$NSP)))
stopifnot(identical(names(seqs$CDS), names(seqs$SP)))

dinuc_levels <- c(
  "AA","AC","AG","AU","CA","CC","CG","CU",
  "GA","GC","GG","GU","UA","UC","UG","UU"
)

dna_to_rna_dinuc <- c(
  AA="AA", AC="AC", AG="AG", AT="AU",
  CA="CA", CC="CC", CG="CG", CT="CU",
  GA="GA", GC="GC", GG="GG", GT="GU",
  TA="UA", TC="UC", TG="UG", TT="UU"
)

calc_rad_one <- function(seq_dna, accession, region) {
  dn <- dinucleotideFrequency(seq_dna, as.prob = FALSE)
  dn <- dn[names(dna_to_rna_dinuc)]
  names(dn) <- dna_to_rna_dinuc[names(dn)]
  dn <- dn[dinuc_levels]

  nf <- alphabetFrequency(seq_dna, baseOnly = TRUE, as.prob = TRUE)
  fA <- nf[["A"]]; fC <- nf[["C"]]; fG <- nf[["G"]]; fT <- nf[["T"]]

  base_freq <- c(A = fA, U = fT, G = fG, C = fC)

  n_pos <- length(seq_dna) - 1
  f_obs <- as.numeric(dn) / n_pos

  first  <- substr(dinuc_levels, 1, 1)
  second <- substr(dinuc_levels, 2, 2)
  f_exp  <- as.numeric(base_freq[first] * base_freq[second])

  tibble(
    Accession    = accession,
    Region       = region,
    Dinucleotide = dinuc_levels,
    Count        = as.integer(dn),
    Freq_obs     = f_obs,
    Freq_exp     = f_exp,
    RAD          = f_obs / f_exp
  )
}

rad_list <- list()
for (region in c("CDS", "NSP", "SP")) {
  s <- seqs[[region]]
  rad_list[[region]] <- purrr::map2_dfr(
    as.list(s), names(s),
    ~ calc_rad_one(.x, .y, region)
  )
}

rad_long <- bind_rows(rad_list)

stopifnot(nrow(rad_long) == 52 * 3 * 16)
stopifnot(!any(is.na(rad_long$RAD)))
stopifnot(!any(is.na(rad_long$Freq_obs)))
stopifnot(!any(is.na(rad_long$Freq_exp)))
stopifnot(all(rad_long$RAD > 0))

for (region in c("CDS", "NSP", "SP")) {
  out <- rad_long %>% filter(Region == region)
  write_csv(out, file.path(out_tab, paste0("RAD_per_sequence_", region, ".csv")))
}

rad_summary <- rad_long %>%
  group_by(Region, Dinucleotide) %>%
  summarise(
    Mean  = mean(RAD),
    SD    = sd(RAD),
    Min   = min(RAD),
    Max   = max(RAD),
    Range = paste0(sprintf("%.3f", min(RAD)), "-", sprintf("%.3f", max(RAD))),
    Status = case_when(
      Mean > 1.23 ~ "Over",
      Mean < 0.78 ~ "Under",
      TRUE        ~ "Neutral"
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Region       = factor(Region, levels = c("CDS", "NSP", "SP")),
    Dinucleotide = factor(Dinucleotide, levels = dinuc_levels)
  ) %>%
  arrange(Region, Dinucleotide)

write_csv(rad_summary, file.path(out_tab, "RAD_summary.csv"))

# ---------- figure ----------
region_colors <- c(CDS = "#C9E7F7", NSP = "#764C7C", SP = "#FFB3E0")

common_theme <- theme_classic(base_size = 10, base_family = "ArialMT") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0, family = "ArialMT"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, color = "black", family = "ArialMT"),
    axis.text.y = element_text(size = 10, color = "black", family = "ArialMT"),
    axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 20), family = "ArialMT"),
    axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 12), family = "ArialMT"),
    legend.position = "none",
    axis.line = element_line(linewidth = 0.75 * 0.352778, colour = "black"),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 16, face = "bold", family = "ArialMT")
  )

p_rad <- ggplot(rad_summary, aes(x = Dinucleotide, y = Mean)) +
  geom_col(aes(fill = Region), width = 0.8, colour = NA) +
  geom_text(aes(label = sprintf("%.3f", Mean), y = Mean + 0.02),
            size = 10 * 0.352778, family = "ArialMT", vjust = 0, hjust = 0.5) +
  geom_hline(yintercept = 1.23, colour = "#E41A1C", linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = 0.78, colour = "#E41A1C", linetype = "dashed", linewidth = 0.5) +
  facet_wrap(~ Region, nrow = 1) +
  scale_fill_manual(values = region_colors) +
  scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.25), expand = c(0, 0)) +
  labs(x = "Dinucleotide", y = "RAD") +
  common_theme

ggsave(file.path(out_fig, "RAD_barplot.pdf"), p_rad,
       width = 18, height = 6, device = "pdf")

cat("Done.\n")
cat("Rows in long table:", nrow(rad_long), "\n")
cat("Files written to:", out_tab, "\n")
cat("Figure written to:", file.path(out_fig, "RAD_barplot.pdf"), "\n")
cat("\nSummary (Region x Dinucleotide):\n")
print(rad_summary %>% select(Region, Dinucleotide, Mean, SD, Status), n = 48)
