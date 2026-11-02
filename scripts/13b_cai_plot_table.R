# ============================================================
# Purpose : Generate Table 4 and the combined CAI figures.
# Input   : Hard-coded data obtained from the CAIcal web server 
#           (Puigbo et al., 2008; http://genomes.urv.cat/CAIcal/).
# Output  : results/tables/CAI/Table4_CAI_web_summary.csv
#           results/figures/CAI_combined_barplots.pdf
# Dependencies: here, tidyverse, patchwork
# Random seed  : no stochastic step
# Run     : Rscript scripts/13b_cai_plot_table.R
# ============================================================

library(here)
library(tidyverse)
library(patchwork)

out_tab <- here("results", "tables", "CAI")
out_fig <- here("results", "figures")
dir.create(out_tab, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

host_levels <- c("Aedes aegypti", "Culex tarsalis", "Haemorhous mexicanus", 
                 "Passer domesticus", "Equus caballus", "Homo sapiens")
region_levels <- c("CDS", "NSP", "SP")

dat <- data.frame(
  Region = factor(rep(region_levels, each = 6), levels = region_levels),
  Host = factor(rep(host_levels, 3), levels = host_levels),
  Average_CAI = c(
    0.742, 0.599, 0.698, 0.732, 0.719, 0.740,
    0.747, 0.605, 0.706, 0.737, 0.723, 0.746,
    0.745, 0.607, 0.705, 0.735, 0.723, 0.744
  ),
  eCAI_threshold = c(
    0.785, 0.658, 0.754, 0.789, 0.772, 0.796,
    0.791, 0.662, 0.766, 0.795, 0.778, 0.803,
    0.788, 0.663, 0.763, 0.792, 0.775, 0.800
  )
)

dat <- dat %>% mutate(CAI_eCAI_Ratio = Average_CAI / eCAI_threshold)

table4 <- dat %>%
  select(Dataset = Region, Host, Observed_CAI = Average_CAI, 
         Expected_eCAI = eCAI_threshold, CAI_eCAI_Ratio) %>%
  mutate(CAI_eCAI_Ratio = round(CAI_eCAI_Ratio, 3))

readr::write_csv(table4, file.path(out_tab, "Table4_CAI_web_summary.csv"))

region_colors <- c(CDS = "#C9E7F7", NSP = "#764C7C", SP = "#FFB3E0")

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

p_cai <- ggplot(dat, aes(x = Host, y = Average_CAI)) +
  geom_bar(aes(fill = Region), stat = "identity", width = 0.8, colour = NA) +
  geom_crossbar(aes(y = eCAI_threshold, ymin = eCAI_threshold, ymax = eCAI_threshold), 
                width = 0.8, colour = "#E41A1C", linetype = "dashed", linewidth = 0.5, fatten = 0) +
  geom_text(aes(label = sprintf("%.3f", Average_CAI), y = Average_CAI + 0.008), 
            size = 10 * 0.352778, family = "ArialMT", vjust = 0, hjust = 0.5) +
  geom_text(aes(label = sprintf("%.3f", eCAI_threshold), y = eCAI_threshold + 0.008), 
            size = 10 * 0.352778, family = "ArialMT", vjust = 0, hjust = 0.5, colour = "#E41A1C") +
  facet_wrap(~ Region, nrow = 1) + labs(title = "A. Codon Adaptation Index (CAI) across hosts", x = NULL, y = "CAI") +
  scale_fill_manual(values = region_colors) +
  scale_x_discrete(labels = function(x) parse(text = paste0("italic('", gsub("_", " ", x), "')"))) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, 0.2), expand = c(0, 0)) +
  common_theme

p_ratio <- ggplot(dat, aes(x = Host, y = CAI_eCAI_Ratio)) +
  geom_bar(aes(fill = Region), stat = "identity", width = 0.8, colour = NA) +
  geom_hline(yintercept = 1.0, linetype = "dashed", colour = "#E41A1C", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", CAI_eCAI_Ratio), y = CAI_eCAI_Ratio + 0.01), 
            size = 10 * 0.352778, family = "ArialMT", vjust = 0, hjust = 0.5) +
  facet_wrap(~ Region, nrow = 1) + labs(title = "B. Normalized CAI (CAI / eCAI Ratio) across hosts", x = NULL, y = "CAI / eCAI Ratio") +
  scale_fill_manual(values = region_colors) +
  scale_x_discrete(labels = function(x) parse(text = paste0("italic('", gsub("_", " ", x), "')"))) +
  scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1.2, 0.2), expand = c(0, 0)) +
  common_theme

combined_plot <- p_cai / p_ratio + plot_layout(guides = "collect") & theme(legend.position = "right")
ggsave(file.path(out_fig, "CAI_combined_barplots.pdf"), combined_plot, width = 18, height = 12, device = "pdf")
