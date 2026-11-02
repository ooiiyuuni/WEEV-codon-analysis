# ============================================================
# 18_coa_permanova.R
# Purpose : PERMANOVA on COA coordinates for CDS, NSP and SP,
#           testing whether host group or geographic region
#           explains codon usage variation.
# Input   : results/tables/COA/COA_coords_{CDS,NSP,SP}.csv
#           results/tables/metadata.csv
# Output  : results/tables/COA/COA_permanova_results.csv
# Dependencies: here, tidyverse, vegan
# Random seed : 12345 (permutations)
# Run     : Rscript scripts/18_coa_permanova.R
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(vegan)
})

set.seed(12345)

coa_dir <- here("results", "tables", "COA")

# ---------- metadata ----------
meta <- read_csv(here("results", "tables", "metadata.csv"), show_col_types = FALSE) %>%
  select(Accession, Host, Country)

host_group <- function(x) {
  x_low <- tolower(x)
  case_when(
    x_low %in% c("culex tarsalis", "mosquito", "culex sp.") ~ "Mosquito",
    x_low %in% c("equine", "horse") ~ "Equine",
    x_low %in% c("homo sapiens") ~ "Human",
    TRUE ~ "Others"
  )
}

region_group <- function(x) {
  prefix <- str_extract(x, "^[^:]+") %>% str_trim()
  case_when(
    prefix %in% c("Uruguay", "Brazil", "Argentina") ~ "South America",
    prefix %in% c("USA", "Canada") ~ "North America",
    TRUE ~ NA_character_
  )
}

meta <- meta %>%
  mutate(
    HostGroup = factor(host_group(Host),
                       levels = c("Mosquito", "Equine", "Human", "Others")),
    GeoRegion = factor(region_group(Country),
                       levels = c("North America", "South America"))
  )
stopifnot(!any(is.na(meta$HostGroup)))
stopifnot(!any(is.na(meta$GeoRegion)))

# ---------- run PERMANOVA ----------
results <- list()

for (r in c("CDS", "NSP", "SP")) {
  d <- read_csv(file.path(coa_dir, paste0("COA_coords_", r, ".csv")),
                show_col_types = FALSE) %>%
    inner_join(meta, by = "Accession")

  stopifnot(nrow(d) == 52)
  stopifnot(!any(is.na(d$Axis1)))
  stopifnot(!any(is.na(d$Axis2)))

  xy <- d %>% select(Axis1, Axis2) %>% as.data.frame()

  for (factor_name in c("HostGroup", "GeoRegion")) {
    f <- d[[factor_name]]
    ad <- adonis2(xy ~ f, method = "euclidean", permutations = 999)

    results[[length(results) + 1]] <- tibble(
      Region = r,
      Factor = factor_name,
      df = ad$Df[1],
      SumSq = ad$SumOfSqs[1],
      R2 = ad$R2[1],
      F = ad$F[1],
      P = ad$`Pr(>F)`[1]
    )
  }
}

res <- bind_rows(results) %>%
  mutate(
    R2 = round(R2, 3),
    F  = round(F, 3),
    P  = signif(P, 3)
  )

write_csv(res, file.path(coa_dir, "COA_permanova_results.csv"))

cat("Done.\n")
cat("PERMANOVA results:\n")
print(res, n = 10)
cat("\nFile written to:", file.path(coa_dir, "COA_permanova_results.csv"), "\n")
