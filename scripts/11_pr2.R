# ============================================================
# 11_pr2.R
# Purpose : PR2 plots (G3/(G3+C3) vs A3/(A3+U3)) for WEEV CDS, NSP
#           and SP, colored by host group.
# Input   : data/WEEV/cds.fasta, nsp.fasta, sp.fasta
#           results/tables/metadata.csv
# Output  : results/figures/PR2_all.pdf
#           results/tables/PR2/PR2_coords_CDS.csv
#           results/tables/PR2/PR2_coords_NSP.csv
#           results/tables/PR2/PR2_coords_SP.csv
# Dependencies: here, Biostrings, tidyverse, patchwork
# Random seed : none (no stochastic step)
# Run     : Rscript scripts/11_pr2.R
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(Biostrings)
  library(tidyverse)
  library(patchwork)
})

fam <- "ArialMT"

tab_dir <- here("results", "tables")
fig_dir <- here("results", "figures")
out_tab <- here("results", "tables", "PR2")

meta <- read.csv(file.path(tab_dir, "metadata.csv"),
                 check.names = FALSE, stringsAsFactors = FALSE) %>%
  select(Accession, Host) %>%
  mutate(meta_row = row_number())

host_group <- function(h) {
  h <- tolower(trimws(h))
  dplyr::case_when(
    h %in% c("equine", "horse") ~ "Equine",
    h %in% c("culex tarsalis", "mosquito", "culex sp.") ~ "Mosquito",
    h %in% c("homo sapiens") ~ "Human",
    TRUE ~ "Others"
  )
}

compute_pr2 <- function(fa_path) {
  dna <- readDNAStringSet(fa_path)
  res <- lapply(seq_along(dna), function(i) {
    s <- toupper(as.character(dna[[i]]))
    n <- nchar(s)
    cod <- substring(s, seq(1, n - 2, by = 3), seq(3, n, by = 3))
    t3 <- substr(cod, 3, 3)
    A3 <- mean(t3 == "A"); U3 <- mean(t3 == "T")
    G3 <- mean(t3 == "G"); C3 <- mean(t3 == "C")
    data.frame(Accession = names(dna)[i],
               x = G3 / (G3 + C3),
               y = A3 / (A3 + U3))
  })
  bind_rows(res) %>%
    mutate(Accession = sub("_(complete_CDS|nsp|sp)$", "", Accession))
}

regions <- c("CDS", "NSP", "SP")
files <- c(CDS = "cds.fasta", NSP = "nsp.fasta", SP = "sp.fasta")

dat_all <- bind_rows(lapply(regions, function(r) {
  d <- compute_pr2(file.path(here("data", "WEEV"), files[[r]]))
  d$Region <- r
  d
})) %>%
  inner_join(meta, by = "Accession") %>%
  mutate(Group = factor(host_group(Host),
                        levels = c("Mosquito", "Equine", "Human", "Others")))

stopifnot(nrow(dat_all) == 156)
for (r in regions) {
  d <- dat_all %>% filter(Region == r)
  stopifnot(nrow(d) == 52)
  stopifnot(!any(is.na(d$Host)), !any(is.na(d$x)), !any(is.na(d$y)))
}
ord <- lapply(regions, function(r) dat_all %>% filter(Region == r) %>% pull(Accession))
stopifnot(identical(ord[[1]], ord[[2]]), identical(ord[[1]], ord[[3]]))

# x axis: left edge fixed at 0.455 (0.46 sits half a 0.01 step from y axis);
# right edge from data with 5% padding, same as before
x_rng <- range(dat_all$x)
x_rng <- c(0.455, x_rng[2] + diff(x_rng) * 0.05)
y_rng <- range(dat_all$y)
y_rng <- y_rng + c(-1, 1) * diff(y_rng) * 0.05

grp_cols <- c(Mosquito = "#64A9E2", Equine = "#DDCE43",
              Human = "#E41A1C", Others = "#681F72")

base_theme <- theme_classic(base_size = 10, base_family = fam) +
  theme(
    plot.title   = element_text(size = 24, face = "bold", hjust = 0.5,
                                family = fam, margin = margin(b = 8)),
    axis.title.x = element_text(size = 24, face = "bold", family = fam,
                                margin = margin(t = 10)),
    axis.title.y = element_text(size = 24, face = "bold", family = fam,
                                margin = margin(r = 20)),
    axis.text    = element_text(size = 10, color = "black", family = fam),
    axis.line    = element_line(color = "black", linewidth = 0.5 * 0.352778),
    axis.ticks   = element_line(color = "black", linewidth = 0.5 * 0.352778),
    legend.title = element_text(size = 10, face = "bold", family = fam),
    legend.text  = element_text(size = 10, family = fam),
    plot.margin  = margin(8, 8, 8, 8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.5 * 0.352778)
  )

make_pr2 <- function(df, x_rng, y_rng, title, show_legend = FALSE) {
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_vline(xintercept = 0.5, color = "grey40", linewidth = 0.5, linetype = "dashed") +
    geom_point(aes(color = Group), shape = 16, size = 2.4, alpha = 0.9) +
    scale_color_manual(values = grp_cols, name = "Host group") +
    scale_x_continuous(limits = x_rng, expand = c(0, 0),
                       breaks = seq(0.46, 0.53, by = 0.01)) +
    scale_y_continuous(limits = y_rng, expand = c(0, 0),
                       breaks = scales::breaks_pretty(n = 5)) +
    coord_fixed(ratio = 1) +
    labs(title = title, x = "G3 / (G3 + C3)", y = "A3 / (A3 + U3)") +
    base_theme
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

plots <- list()
for (r in regions) {
  plots[[r]] <- make_pr2(dat_all %>% filter(Region == r), x_rng, y_rng, r,
                         show_legend = (r == "CDS"))
}
combined <- (plots[["CDS"]] | plots[["NSP"]] | plots[["SP"]]) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "PR2_all.pdf"),
       combined, width = 18, height = 6, device = "pdf")

for (r in regions) {
  d <- dat_all %>% filter(Region == r) %>%
    select(Accession, Host, Group, x, y, meta_row) %>%
    arrange(meta_row) %>% select(-meta_row)

  mean_row <- data.frame(Accession = "Mean", Host = "", Group = "",
                         x = mean(d$x), y = mean(d$y))
  sd_row   <- data.frame(Accession = "SD", Host = "", Group = "",
                         x = sd(d$x), y = sd(d$y))

  out <- bind_rows(d, mean_row, sd_row) %>%
    mutate(x = sprintf("%.4f", x),
           y = sprintf("%.4f", y))

  write.csv(out, file.path(out_tab, paste0("PR2_coords_", r, ".csv")),
            row.names = FALSE, quote = FALSE)
}

cat("Done.\n")
cat("Figure:", file.path(fig_dir, "PR2_all.pdf"), "\n")
cat("Tables:", out_tab, "\n")
