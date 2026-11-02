# ============================================================
# 12_neutrality.R
# Purpose : Neutrality plots (GC12 vs GC3) for WEEV CDS, NSP and SP,
#           colored by host group; linear regression per region with
#           slope, intercept, R2 and P value.
# Input   : data/WEEV/cds.fasta, nsp.fasta, sp.fasta
#           results/tables/metadata.csv
# Output  : results/figures/Neutrality_all.pdf
#           results/tables/Neutrality/Neutrality_per_sequence_{CDS,NSP,SP}.csv
#           results/tables/Neutrality/Neutrality_regression.csv
# Dependencies: here, Biostrings, tidyverse, patchwork
# Random seed : none (no stochastic step)
# Run     : Rscript scripts/12_neutrality.R
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
out_tab <- here("results", "tables", "Neutrality")

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

compute_gc <- function(fa_path) {
  dna <- readDNAStringSet(fa_path)
  res <- lapply(seq_along(dna), function(i) {
    s <- toupper(as.character(dna[[i]]))
    n <- nchar(s)
    cod <- substring(s, seq(1, n - 2, by = 3), seq(3, n, by = 3))
    GC1 <- mean(substr(cod, 1, 1) %in% c("G","C"))
    GC2 <- mean(substr(cod, 2, 2) %in% c("G","C"))
    GC3 <- mean(substr(cod, 3, 3) %in% c("G","C"))
    data.frame(Accession = names(dna)[i], GC12 = (GC1 + GC2) / 2, GC3 = GC3)
  })
  bind_rows(res) %>%
    mutate(Accession = sub("_(complete_CDS|nsp|sp)$", "", Accession))
}

regions <- c("CDS", "NSP", "SP")
files <- c(CDS = "cds.fasta", NSP = "nsp.fasta", SP = "sp.fasta")

dat_all <- bind_rows(lapply(regions, function(r) {
  d <- compute_gc(file.path(here("data", "WEEV"), files[[r]]))
  d$Region <- r
  d
})) %>%
  inner_join(meta, by = "Accession") %>%
  mutate(Group = factor(host_group(Host),
                        levels = c("Mosquito", "Equine", "Human", "Others")))

# ---- sanity checks ----
stopifnot(nrow(dat_all) == 156)
for (r in regions) {
  d <- dat_all %>% filter(Region == r)
  stopifnot(nrow(d) == 52)
  stopifnot(!any(is.na(d$Host)), !any(is.na(d$GC3)), !any(is.na(d$GC12)))
}
ord <- lapply(regions, function(r) dat_all %>% filter(Region == r) %>% pull(Accession))
stopifnot(identical(ord[[1]], ord[[2]]), identical(ord[[1]], ord[[3]]))

# ---- regression per region ----
fit_one <- function(d) {
  fit <- lm(GC12 ~ GC3, data = d)
  s   <- summary(fit)
  data.frame(slope = unname(coef(fit)["GC3"]),
             intercept = unname(coef(fit)["(Intercept)"]),
             R2 = s$r.squared,
             P_value = s$coefficients["GC3", "Pr(>|t|)"],
             n = nrow(d))
}
reg_tab <- bind_rows(lapply(regions, function(r) {
  fit_one(dat_all %>% filter(Region == r)) %>% mutate(region = r, .before = 1)
}))

# ---- shared axis ranges across regions, padded 20% ----
x_rng <- range(dat_all$GC3)
x_rng <- x_rng + c(-1, 1) * diff(x_rng) * 0.20
y_rng <- range(dat_all$GC12)
y_rng <- y_rng + c(-1, 1) * diff(y_rng) * 0.20

grp_cols <- c(Mosquito = "#64A9E2", Equine = "#DDCE43",
              Human = "#E41A1C", Others = "#681F72")

base_theme <- theme_classic(base_size = 10, base_family = fam) +
  theme(
    plot.title   = element_text(size = 16, face = "bold", hjust = 0.5,
                                family = fam, margin = margin(b = 8)),
    axis.title.x = element_text(size = 16, face = "bold", family = fam,
                                margin = margin(t = 10)),
    axis.title.y = element_text(size = 16, face = "bold", family = fam,
                                margin = margin(r = 20)),
    axis.text    = element_text(size = 10, color = "black", family = fam),
    axis.line    = element_line(color = "black", linewidth = 0.75 * 0.352778),
    axis.ticks   = element_line(color = "black", linewidth = 0.75 * 0.352778),
    legend.title = element_text(size = 10, face = "bold", family = fam),
    legend.text  = element_text(size = 10, family = fam),
    plot.margin  = margin(8, 8, 8, 8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3)
  )

make_plot <- function(d, reg, x_rng, y_rng, title, show_legend = FALSE) {
  lbl <- sprintf("y = %.4fx + %.4f\nR2 = %.3f\nP = %.3f\nn = %d",
                 reg$slope, reg$intercept, reg$R2, reg$P_value, reg$n)
  p <- ggplot(d, aes(x = GC3, y = GC12)) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                color = "black", linewidth = 0.8) +
    geom_point(aes(color = Group), shape = 16, size = 2.4, alpha = 0.9) +
    scale_color_manual(values = grp_cols, name = "Host group") +
    scale_x_continuous(limits = x_rng, expand = c(0, 0),
                       breaks = scales::breaks_pretty(n = 5)) +
    scale_y_continuous(limits = y_rng, expand = c(0, 0),
                       breaks = scales::breaks_pretty(n = 5)) +
    coord_fixed(ratio = (x_rng[2] - x_rng[1]) / (y_rng[2] - y_rng[1])) +
    annotate("text", x = x_rng[1] + 0.02 * diff(x_rng),
             y = y_rng[2] - 0.02 * diff(y_rng),
             label = lbl, hjust = 0, vjust = 1, size = 3, family = fam) +
    labs(title = paste0(title, " neutrality plot"), x = "GC3", y = "GC12") +
    base_theme
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

plots <- list()
for (r in regions) {
  d <- dat_all %>% filter(Region == r)
  reg <- reg_tab %>% filter(region == r)
  plots[[r]] <- make_plot(d, reg, x_rng, y_rng, r, show_legend = (r == "CDS"))
}
combined <- (plots[["CDS"]] | plots[["NSP"]] | plots[["SP"]]) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "Neutrality_all.pdf"),
       combined, width = 18, height = 5, device = "pdf")

# ---- per-sequence tables with Mean and SD rows ----
for (r in regions) {
  d <- dat_all %>% filter(Region == r) %>%
    select(Accession, Host, Group, GC3, GC12, meta_row) %>%
    arrange(meta_row) %>% select(-meta_row)
  mean_row <- data.frame(Accession = "Mean", Host = "", Group = "",
                         GC3 = mean(d$GC3), GC12 = mean(d$GC12))
  sd_row   <- data.frame(Accession = "SD", Host = "", Group = "",
                         GC3 = sd(d$GC3), GC12 = sd(d$GC12))
  out <- bind_rows(d, mean_row, sd_row) %>%
    mutate(GC3 = sprintf("%.4f", GC3), GC12 = sprintf("%.4f", GC12))
  write.csv(out, file.path(out_tab, paste0("Neutrality_per_sequence_", r, ".csv")),
            row.names = FALSE, quote = FALSE)
}

write.csv(reg_tab, file.path(out_tab, "Neutrality_regression.csv"), row.names = FALSE)

cat("Done.\n")
cat("Figure:", file.path(fig_dir, "Neutrality_all.pdf"), "\n")
cat("Tables:", out_tab, "\n")
print(reg_tab, row.names = FALSE)
