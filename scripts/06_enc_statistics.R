# ============================================================
# 06_enc_statistics.R
# Purpose : Compare ENc among the complete CDS, NSP, and SP datasets
# Method  : Kruskal-Wallis test + Dunn's post hoc test (BH correction)
# Input   : results/tables/codonw/CodonUsageIndices_{CDS,NSP,SP}.csv
# Output  : results/tables/ENc_KruskalWallis_summary.csv
#           results/tables/ENc_Dunn_posthoc.csv
# Depends : here, tidyverse, rstatix (v1.1.0), R (v4.6.1)
# Seed    : none
# Run     : Rscript scripts/06_enc_statistics.R
# ============================================================

library(here)
library(tidyverse)
library(rstatix)

# Read one dataset and keep only ENc, dropping summary rows.
read_enc <- function(dataset) {
  read_csv(
    here("results", "tables", "codonw",
         paste0("CodonUsageIndices_", dataset, ".csv")),
    show_col_types = FALSE
  ) %>%
    filter(!str_detect(Accession, "_Mean$|_SD$")) %>%
    select(Accession, ENc = Nc) %>%
    mutate(dataset = dataset)
}

enc_all <- bind_rows(
  read_enc("CDS"),
  read_enc("NSP"),
  read_enc("SP")
)

# Expect 52 sequences per dataset, 156 total.
stopifnot(nrow(enc_all) == 156)
stopifnot(all(table(enc_all$dataset) == 52))

# Global test
kw <- kruskal.test(ENc ~ dataset, data = enc_all)
print(kw)

# Pairwise post hoc with BH correction
dunn <- dunn_test(enc_all, ENc ~ dataset, p.adjust.method = "BH")
print(dunn)

# Write outputs
kw_summary <- tibble(
  test        = "Kruskal-Wallis",
  chi_squared = unname(kw$statistic),
  df          = unname(kw$parameter),
  p_value     = kw$p.value
)

write_csv(kw_summary, here("results", "tables", "ENc_KruskalWallis_summary.csv"))
write_csv(dunn,       here("results", "tables", "ENc_Dunn_posthoc.csv"))