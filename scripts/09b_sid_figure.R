# =============================================================================
# 09b_sid_figure.R
#
# Purpose : Summarise mean SiD by host and genomic region, write the summary
#           table, and generate a grouped bar chart of mean SiD between WEEV
#           sequences and six host reference profiles. Three genomic regions
#           (CDS, NSP, SP) are shown side by side for each host, with bars
#           flush within each host group. Bar-top labels are horizontal and
#           placed at a fixed absolute offset above every bar. Output is SVG
#           so that every bar and text label remains an individually
#           selectable object in Adobe Illustrator.
#
# Input   : results/tables/SiD/SiD_long.csv
#           (columns: Accession, Region, Host, R, SiD)
#
# Output  : results/tables/SiD/SiD_summary.csv  (18 rows)
#           results/figures/SiD_barplot.svg
#
# Dependencies: here, tidyverse, svglite
#
# Random seed : set.seed(12345) (no stochastic step; project convention)
#
# Run     : Rscript scripts/09b_sid_figure.R
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(svglite)
})

set.seed(12345)

# ---------------------------------------------------------------------------
# 1. Read per-sequence SiD
# ---------------------------------------------------------------------------

sid <- read_csv(here("results/tables/SiD/SiD_long.csv"), show_col_types = FALSE)

# ---------------------------------------------------------------------------
# 2. Fixed host order and region order
# ---------------------------------------------------------------------------

host_order <- c(
  "Aedes aegypti",
  "Culex tarsalis",
  "Haemorhous mexicanus",
  "Passer domesticus",
  "Equus caballus",
  "Homo sapiens"
)

region_order <- c("CDS", "NSP", "SP")

# ---------------------------------------------------------------------------
# 3. Summarise mean SiD by Host and Region, write summary table
# ---------------------------------------------------------------------------

sid_summary <- sid %>%
  mutate(
    Host   = factor(Host,   levels = host_order),
    Region = factor(Region, levels = region_order)
  ) %>%
  group_by(Host, Region) %>%
  summarise(mean_SiD = mean(SiD), .groups = "drop")

stopifnot(nrow(sid_summary) == 18)

dir.create(here("results/tables/SiD"), recursive = TRUE, showWarnings = FALSE)
write_csv(sid_summary, here("results/tables/SiD/SiD_summary.csv"))

# ---------------------------------------------------------------------------
# 4. Plot
# ---------------------------------------------------------------------------

region_colours <- c(
  "CDS" = "#C9E7F7",
  "NSP" = "#764C7C",
  "SP"  = "#FFB3E0"
)

bar_width   <- 0.95
dodge_width <- 0.95

# Fixed absolute offset of the label above each bar top (3% of y range)
label_nudge <- max(sid_summary$mean_SiD) * 0.03

p <- ggplot(sid_summary, aes(x = Host, y = mean_SiD, fill = Region)) +
  geom_col(
    position  = position_dodge(width = dodge_width),
    width     = bar_width,
    colour    = NA,
    linewidth = 0
  ) +
  geom_text(
    aes(
      label = formatC(mean_SiD, format = "f", digits = 4),
      group = Region
    ),
    position  = position_dodge(width = dodge_width),
    vjust     = 0,
    nudge_y   = label_nudge,
    size      = 7 / .pt,
    family    = "Arial"
  ) +
  scale_fill_manual(values = region_colours, name = "Region") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  coord_cartesian(clip = "off") +
  labs(
    title = "SiD between WEEV and host species",
    x     = "Host species",
    y     = "SiD"
  ) +
  theme_classic(base_size = 10, base_family = "Arial") +
  theme(
    plot.title        = element_text(size = 16, face = "bold", hjust = 0.5,
                                     margin = margin(b = 12)),
    axis.title.x      = element_text(size = 16, face = "bold", margin = margin(t = 12)),
    axis.title.y      = element_text(size = 16, face = "bold", angle = 90,
                                     margin = margin(r = 30)),
    axis.text.x       = element_text(size = 10, face = "italic", colour = "black",
                                     margin = margin(t = 4)),
    axis.text.y       = element_text(size = 10, colour = "black"),
    axis.line         = element_line(colour = "black", linewidth = 0.75 * 0.352778),
    axis.ticks        = element_line(colour = "black", linewidth = 0.75 * 0.352778),
    legend.title      = element_text(size = 10),
    legend.text       = element_text(size = 10),
    legend.key        = element_rect(fill = NA, colour = NA),
    legend.position   = "right",
    panel.grid        = element_blank(),
    plot.background   = element_blank(),
    panel.background  = element_blank(),
    plot.margin       = margin(10, 10, 10, 10)
  )

# ---------------------------------------------------------------------------
# 5. Save as editable SVG
# ---------------------------------------------------------------------------

svglite(
  filename = here("results/figures/SiD_barplot.svg"),
  width    = 9,
  height   = 5.5
)
print(p)
dev.off()

message("SiD_summary.csv and SiD_barplot.svg written")