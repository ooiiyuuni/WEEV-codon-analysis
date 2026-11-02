# ============================================================
# Purpose : Generate the RCDI table and figure.
# Input   : Hard-coded data obtained from the RCDI/eRCDI web server 
#           (Puigbo et al., 2010; http://genomes.urv.cat/RCDI/).
# Output  : results/tables/RCDI/Table_RCDI_web_summary.csv
#           results/figures/RCDI_barplot.pdf
# Dependencies: here, tidyverse
# Random seed  : no stochastic step
# Run     : Rscript scripts/14_rcdi_statistics.R
# ============================================================

library(here)
library(tidyverse)

out_tab <- here("results", "tables", "RCDI")
out_fig <- here("results", "figures")
dir.create(out_tab, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

host_levels <- c("Aedes aegypti", "Culex tarsalis", "Haemorhous mexicanus", 
                 "Passer domesticus", "Equus caballus", "Homo sapiens")
region_levels <- c("NSP", "SP")

dat <- data.frame(
  Region = factor(rep(region_levels, each = 6), levels = region_levels),
  Host = factor(rep(host_levels, 2), levels = host_levels),
  Average_RCDI = c(
    1.228, 1.442, 1.340, 1.319, 1.293, 1.310,
    1.230, 1.419, 1.325, 1.314, 1.291, 1.316
  ),
  eRCDI_threshold = c(
    1.312, 1.609, 1.474, 1.449, 1.405, 1.436,
    1.323, 1.577, 1.447, 1.445, 1.409, 1.451
  )
)

readr::write_csv(dat, file.path(out_tab, "Table_RCDI_web_summary.csv"))

region_colors <- c(NSP = "#764C7C", SP = "#FFB3E0")

common_theme <- theme_classic(base_size = 10, base_family = "ArialMT") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0, margin = margin(b = 10), family = "ArialMT"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 10, color = "black", family = "ArialMT"),
    axis.text.y = element_text(size = 10, color = "black", family = "ArialMT"),
    axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 12), family = "ArialMT"),
    legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 10, family = "ArialMT"),
    axis.line = element_line(linewidth = 0.75 * 0.352778, colour = "black"),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3), panel.grid.minor = element_blank(),
    strip.background = element_blank(), strip.text = element_text(size = 10, face = "bold", family = "ArialMT")
  )

p_rcdi <- ggplot(dat, aes(x = Host, y = Average_RCDI)) +
  geom_bar(aes(fill = Region), stat = "identity", width = 0.8, colour = NA) +
  geom_crossbar(aes(y = eRCDI_threshold, ymin = eRCDI_threshold, ymax = eRCDI_threshold), 
                width = 0.8, colour = "#E41A1C", linetype = "dashed", linewidth = 0.5, fatten = 0) +
  geom_text(aes(label = sprintf("%.3f", Average_RCDI), y = Average_RCDI + 0.02), 
            size = 10 * 0.352778, family = "ArialMT", vjust = 0, hjust = 0.5) +
  geom_text(aes(label = sprintf("%.3f", eRCDI_threshold), y = eRCDI_threshold + 0.02), 
            size = 10 * 0.352778, family = "ArialMT", vjust = 0, hjust = 0.5, colour = "#E41A1C") +
  facet_wrap(~ Region, nrow = 1) +
  labs(title = "Relative Codon Deoptimization Index (RCDI) across hosts", x = NULL, y = "RCDI") +
  scale_fill_manual(values = region_colors) +
  scale_x_discrete(labels = function(x) parse(text = paste0("italic('", gsub("_", " ", x), "')"))) +
  scale_y_continuous(limits = c(0, 2.0), breaks = seq(0, 2.0, 0.5), expand = c(0, 0)) +
  common_theme

ggsave(file.path(out_fig, "RCDI_barplot.pdf"), p_rcdi, width = 18, height = 7, device = "pdf")
