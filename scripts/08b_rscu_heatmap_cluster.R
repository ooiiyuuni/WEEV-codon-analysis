# ============================================================
# 08b_rscu_heatmap_cluster.R
# Purpose: RSCU heatmap with codon clustering.
#          - Hierarchical clustering (ward.D2) of 59 codons based on
#            RSCU vectors across nine datasets
#          - No amino acid grouping bars
#          - Each codon labelled with its own amino acid (vertical text)
#          - Bottom row: five-class colour strip, one tile per codon
#          - Vertical layout matches 08_rscu_heatmap.R
#            (codon y = 0, amino acid y = -1.8, class strip y = -2.6)
#          - Discretized colour scale (blue-white-red, RdBu reversed)
#          - Host names italic, ordered phylogenetically
#          - Arial 10 pt; axis/ticks/border 0.75 pt black
# Input:
#   results/tables/RSCU/RSCU.csv
# Output:
#   results/figures/RSCU_heatmap_cluster.pdf
# Dependencies: here, dplyr, tidyr, readr, ggplot2, ggnewscale
# Random seed: 12345 (hierarchical clustering is deterministic;
#              seed kept per project convention)
# Run: Rscript scripts/08b_rscu_heatmap_cluster.R
# ============================================================

set.seed(12345)  # random seed: 12345

suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr)
  library(readr); library(ggplot2); library(ggnewscale)
})

PT <- 0.352778  # 1 pt = 0.352778 mm

# ---------- 1. Read RSCU table ----------
rscu <- read_csv(here("results","tables","RSCU","RSCU.csv"),
                 show_col_types = FALSE,
                 col_types = cols(.default = col_character()))
num_cols <- setdiff(names(rscu), c("AA","Codon"))
rscu[num_cols] <- lapply(rscu[num_cols], as.numeric)

rscu59 <- rscu %>% filter(!Codon %in% c("AUG","UGG"))

# ---------- 2. Datasets ----------
data_cols <- c("Complete","NSP","SP",
               "Aedes aegypti","Culex tarsalis",
               "Haemorhous mexicanus","Passer domesticus",
               "Equus caballus","Homo sapiens")

# ---------- 3. Hierarchical clustering of codons ----------
mat <- rscu59 %>% select(all_of(data_cols)) %>% as.matrix()
rownames(mat) <- rscu59$Codon
hc <- hclust(dist(mat), method = "ward.D2")
codon_order <- rownames(mat)[hc$order]

# ---------- 4. Five chemical classes ----------
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
class_colors <- c(
  "Nonpolar aliphatic" = "#E6AB02",
  "Aromatic"           = "#7570B3",
  "Polar uncharged"    = "#1B9E77",
  "Positively charged" = "#D95F02",
  "Negatively charged" = "#4D4D4D"
)

rscu59 <- rscu59 %>%
  mutate(Codon = factor(Codon, levels = codon_order),
         CodonIdx = as.integer(Codon),
         Class = factor(aa_class[AA], levels = class_order)) %>%
  arrange(CodonIdx)

# ---------- 5. Long format ----------
long <- rscu59 %>%
  select(Codon, CodonIdx, AA, all_of(data_cols)) %>%
  pivot_longer(cols = all_of(data_cols), names_to = "Dataset", values_to = "RSCU")
long$Dataset <- factor(long$Dataset, levels = rev(data_cols))

# ---------- 6. Discretized colour scale ----------
long <- long %>%
  mutate(RSCU_cut = cut(RSCU,
                        breaks = c(-Inf, 0.6, 0.8, 1.0, 1.2, 1.6, Inf),
                        labels = c("<0.6","0.6-0.8","0.8-1.0","1.0-1.2","1.2-1.6",">1.6"),
                        right = FALSE))
fill_colors <- c("<0.6"="#053061","0.6-0.8"="#92C5DE","0.8-1.0"="#F7F7F7",
                 "1.0-1.2"="#FDDBC7","1.2-1.6"="#D6604D",">1.6"="#67001F")

# ---------- 7. Y-axis labels ----------
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

# ---------- 8. Plot ----------
n_ds <- length(data_cols)
y_codon <- 0
y_aa    <- -1.8
y_class <- -2.6

p <- ggplot(long, aes(x = CodonIdx, y = Dataset, fill = RSCU_cut)) +
  geom_tile(color = NA) +
  scale_fill_manual(values = fill_colors, name = "RSCU", drop = FALSE) +
  geom_text(data = rscu59, aes(x = CodonIdx, y = y_codon, label = Codon),
            inherit.aes = FALSE, size = 10/2.845, angle = 90, family = "Arial",
            color = "black", vjust = 0.5, hjust = 0.5) +
  geom_text(data = rscu59, aes(x = CodonIdx, y = y_aa, label = AA),
            inherit.aes = FALSE, size = 10/2.845, angle = 90, family = "Arial",
            color = "black", vjust = 0.5, hjust = 0.5) +
  new_scale_fill() +
  geom_tile(data = rscu59,
            aes(x = CodonIdx, y = y_class, fill = Class),
            inherit.aes = FALSE, height = 0.8) +
  scale_fill_manual(values = class_colors, name = "Amino acid class") +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  scale_y_discrete(expand = c(0, 0), labels = dataset_labels) +
  coord_cartesian(ylim = c(y_class - 1.6, n_ds + 0.5), clip = "off") +
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

# ---------- 9. Export ----------
out_dir <- here("results","figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(here(out_dir, "RSCU_heatmap_cluster.pdf"),
       p, width = 24, height = 7, device = cairo_pdf)

cat("\n===== RSCU clustered heatmap completed =====\n")
cat("Output: results/figures/RSCU_heatmap_cluster.pdf\n")