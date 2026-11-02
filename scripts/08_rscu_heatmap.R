# ============================================================
# 08_rscu_heatmap.R
# Purpose: RSCU heatmap grouped by amino acid chemical class.
#          Codon order within each synonymous family and amino acid
#          order within each chemical class follow RSCU.csv.
#          Five mutually exclusive classes based on amino acid chemical
#          properties. Discretized colour scale (blue-white-red).
#          Host names italic and ordered phylogenetically.
# Input:
#   results/tables/RSCU/RSCU.csv
# Output:
#   results/figures/RSCU_heatmap.pdf
# Dependencies: here, dplyr, tidyr, readr, ggplot2
# Random seed: 12345 (no stochastic step; kept per project convention)
# Run: Rscript scripts/08_rscu_heatmap.R
# ============================================================

set.seed(12345)  # random seed: 12345

suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr)
  library(readr); library(ggplot2)
})

PT <- 0.352778  # 1 pt = 0.352778 mm

# ---------- 1. Read RSCU table ----------
rscu <- read_csv(here("results","tables","RSCU","RSCU.csv"),
                 show_col_types = FALSE,
                 col_types = cols(.default = col_character()))
num_cols <- setdiff(names(rscu), c("AA","Codon"))
rscu[num_cols] <- lapply(rscu[num_cols], as.numeric)

rscu59 <- rscu %>% filter(!Codon %in% c("AUG","UGG"))

# ---------- 2. Five chemical classes ----------
aa_class <- c(
  "Gly"="Nonpolar aliphatic","Ala"="Nonpolar aliphatic","Val"="Nonpolar aliphatic",
  "Leu"="Nonpolar aliphatic","Ile"="Nonpolar aliphatic","Pro"="Nonpolar aliphatic",
  "Phe"="Aromatic","Tyr"="Aromatic",
  "Ser"="Polar uncharged","Thr"="Polar uncharged","Cys"="Polar uncharged",
  "Asn"="Polar uncharged","Gln"="Polar uncharged",
  "Lys"="Positively charged","Arg"="Positively charged","His"="Positively charged",
  "Asp"="Negatively charged","Glu"="Negatively charged"
)
class_order <- c("Nonpolar aliphatic","Aromatic","Polar uncharged",
                 "Positively charged","Negatively charged")
# Amino acid order within each class follows RSCU.csv
aa_order <- c("Leu","Pro","Ile","Val","Ala","Gly",
              "Phe","Tyr",
              "Ser","Cys","Gln","Thr","Asn",
              "His","Arg","Lys",
              "Asp","Glu")

rscu59 <- rscu59 %>%
  mutate(Class = factor(aa_class[AA], levels = class_order),
         AA = factor(AA, levels = aa_order),
         Codon_orig = factor(Codon, levels = Codon)) %>%
  arrange(Class, AA, Codon_orig) %>%
  mutate(CodonIdx = row_number(),
         Codon = factor(Codon, levels = Codon))

# ---------- 3. Long format ----------
data_cols <- c("Complete","NSP","SP",
               "Aedes aegypti","Culex tarsalis",
               "Haemorhous mexicanus","Passer domesticus",
               "Equus caballus","Homo sapiens")
dataset_levels <- rev(data_cols)

long <- rscu59 %>%
  select(Codon, CodonIdx, AA, Class, all_of(data_cols)) %>%
  pivot_longer(cols = all_of(data_cols), names_to = "Dataset", values_to = "RSCU")
long$Dataset <- factor(long$Dataset, levels = dataset_levels)

# ---------- 4. Discretized colour scale ----------
long <- long %>%
  mutate(RSCU_cut = cut(RSCU,
                        breaks = c(-Inf, 0.6, 0.8, 1.0, 1.2, 1.6, Inf),
                        labels = c("<0.6","0.6-0.8","0.8-1.0","1.0-1.2","1.2-1.6",">1.6"),
                        right = FALSE))
fill_colors <- c("<0.6"="#053061","0.6-0.8"="#92C5DE","0.8-1.0"="#F7F7F7",
                 "1.0-1.2"="#FDDBC7","1.2-1.6"="#D6604D",">1.6"="#67001F")

# ---------- 5. Group bars ----------
gap <- 0.2
aa_bars <- rscu59 %>% group_by(Class, AA) %>%
  summarise(start = min(CodonIdx) - 0.5 + gap,
            end   = max(CodonIdx) + 0.5 - gap,
            mid   = mean(CodonIdx), .groups = "drop")
class_bars <- rscu59 %>% group_by(Class) %>%
  summarise(start = min(CodonIdx) - 0.5 + gap,
            end   = max(CodonIdx) + 0.5 - gap,
            mid   = mean(CodonIdx), .groups = "drop")

# ---------- 6. Y-axis labels ----------
dataset_labels <- c(
  "Complete"             = "CDS",
  "NSP"                  = "NSP",
  "SP"                   = "SP",
  "Aedes aegypti"        = expression(italic("Aedes aegypti")),
  "Culex tarsalis"       = expression(italic("Culex tarsalis")),
  "Haemorhous mexicanus" = expression(italic("Haemorhous mexicanus")),
  "Passer domesticus"    = expression(italic("Passer domesticus")),
  "Equus caballus"       = expression(italic("Equus caballus")),
  "Homo sapiens"         = expression(italic("Homo sapiens"))
)

# ---------- 7. Plot ----------
n_ds <- length(data_cols)
y_codon   <- 0
y_aa_bar  <- -1.0
y_aa_text <- -1.8
y_cls_bar <- -2.6
y_cls_text<- -3.4

p <- ggplot(long, aes(x = CodonIdx, y = Dataset, fill = RSCU_cut)) +
  geom_tile(color = NA) +
  scale_fill_manual(values = fill_colors, name = "RSCU", drop = FALSE) +
  geom_text(data = rscu59, aes(x = CodonIdx, y = y_codon, label = Codon),
            inherit.aes = FALSE, size = 10/2.845, angle = 90, family = "Arial",
            color = "black", vjust = 0.5, hjust = 0.5) +
  geom_segment(data = aa_bars,
               aes(x = start, xend = end, y = y_aa_bar, yend = y_aa_bar),
               inherit.aes = FALSE, linewidth = 0.75 * PT, color = "black") +
  geom_segment(data = class_bars,
               aes(x = start, xend = end, y = y_cls_bar, yend = y_cls_bar),
               inherit.aes = FALSE, linewidth = 0.75 * PT, color = "black") +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  scale_y_discrete(expand = c(0, 0), labels = dataset_labels) +
  coord_cartesian(ylim = c(y_cls_text - 0.8, n_ds + 0.5), clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 10, base_family = "Arial") +
  theme(
    text          = element_text(family = "Arial", size = 10, color = "black"),
    axis.text.y   = element_text(family = "Arial", size = 10, color = "black"),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.ticks.y  = element_line(color = "black", linewidth = 0.75 * PT),
    axis.line     = element_line(color = "black", linewidth = 0.75 * PT),
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 0.75 * PT),
    panel.grid    = element_blank(),
    legend.key    = element_rect(color = "black", linewidth = 0.75 * PT)
  )

for (i in seq_len(nrow(aa_bars))) {
  p <- p + geom_text(data = aa_bars[i, ],
                     aes(x = mid, y = y_aa_text, label = AA),
                     inherit.aes = FALSE, size = 10/2.845,
                     fontface = "plain", family = "Arial", color = "black")
}
for (i in seq_len(nrow(class_bars))) {
  p <- p + geom_text(data = class_bars[i, ],
                     aes(x = mid, y = y_cls_text, label = Class),
                     inherit.aes = FALSE, size = 10/2.845,
                     fontface = "bold", family = "Arial", color = "black")
}

# ---------- 8. Export ----------
out_dir <- here("results","figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(here(out_dir, "RSCU_heatmap.pdf"),
       p, width = 24, height = 7, device = cairo_pdf)

cat("\n===== RSCU heatmap completed =====\n")
cat("Output: results/figures/RSCU_heatmap.pdf\n")