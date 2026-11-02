# ============================================================
# 17_coa_clustering.R
# Purpose : Test whether COA ordination space shows discrete
#           clustering, using k-means (k = 2, 3, 4) and silhouette
#           widths for CDS, NSP and SP. Also tabulate host-group
#           and geographic composition of each cluster and compare
#           cluster membership across the three regions.
# Input   : results/tables/COA/COA_coords_{CDS,NSP,SP}.csv
#           results/tables/metadata.csv
# Output  : results/tables/COA/COA_cluster_silhouette.csv
#           results/tables/COA/COA_cluster_membership.csv
#           results/tables/COA/COA_cluster_host_geo_counts.csv
#           results/tables/COA/COA_cluster_cross_region.csv
# Dependencies: here, tidyverse, cluster
# Random seed : 12345 (k-means starts)
# Run     : Rscript scripts/17_coa_clustering.R
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(cluster)
})

set.seed(12345)

coa_dir <- here("results", "tables", "COA")

# ---------- metadata with host / geo groups ----------
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

# ---------- read coordinates ----------
coords <- lapply(c(CDS = "CDS", NSP = "NSP", SP = "SP"), function(r) {
  read_csv(file.path(coa_dir, paste0("COA_coords_", r, ".csv")), show_col_types = FALSE)
})

for (nm in names(coords)) {
  stopifnot(nrow(coords[[nm]]) == 52)
  stopifnot(!any(is.na(coords[[nm]]$Axis1)))
  stopifnot(!any(is.na(coords[[nm]]$Axis2)))
}

# ---------- k-means and silhouette ----------
ks <- 2:4

sil_tab <- list()
membership <- list()

for (r in names(coords)) {
  d <- coords[[r]]
  xy <- as.matrix(d[, c("Axis1", "Axis2")])
  rownames(xy) <- d$Accession

  for (k in ks) {
    km <- kmeans(xy, centers = k, nstart = 50)
    sil <- silhouette(km$cluster, dist(xy))
    mean_sil <- mean(sil[, "sil_width"])

    sil_tab[[length(sil_tab) + 1]] <- tibble(
      Region = r,
      k = k,
      Mean_silhouette = mean_sil
    )

    membership[[length(membership) + 1]] <- tibble(
      Region = r,
      k = k,
      Accession = rownames(xy),
      Cluster = km$cluster,
      Sil_width = sil[, "sil_width"]
    )
  }
}

sil_tab <- bind_rows(sil_tab)
membership <- bind_rows(membership)

# ---------- write silhouette table ----------
sil_out <- sil_tab %>%
  mutate(Mean_silhouette = round(Mean_silhouette, 3)) %>%
  pivot_wider(names_from = Region, values_from = Mean_silhouette) %>%
  select(k, CDS, NSP, SP)

write_csv(sil_out, file.path(coa_dir, "COA_cluster_silhouette.csv"))

# ---------- attach metadata ----------
membership <- membership %>%
  left_join(meta, by = "Accession")
stopifnot(nrow(membership) == 52 * 3 * 3)
stopifnot(!any(is.na(membership$HostGroup)))

write_csv(membership, file.path(coa_dir, "COA_cluster_membership.csv"))

# ---------- host / geo composition per cluster (k = 3) ----------
counts <- membership %>%
  filter(k == 3) %>%
  group_by(Region, k, Cluster) %>%
  summarise(
    n = n(),
    Mosquito = sum(HostGroup == "Mosquito"),
    Equine   = sum(HostGroup == "Equine"),
    Human    = sum(HostGroup == "Human"),
    Others   = sum(HostGroup == "Others"),
    North_America = sum(GeoRegion == "North America"),
    South_America = sum(GeoRegion == "South America"),
    .groups = "drop"
  )

write_csv(counts, file.path(coa_dir, "COA_cluster_host_geo_counts.csv"))

# ---------- cross-region cluster membership (k = 3) ----------
cross <- membership %>%
  filter(k == 3) %>%
  select(Region, Accession, Cluster) %>%
  pivot_wider(names_from = Region, values_from = Cluster, names_prefix = "Cluster_")

write_csv(cross, file.path(coa_dir, "COA_cluster_cross_region.csv"))

# ---------- console report ----------
cat("Done.\n")
cat("Silhouette widths (mean, k = 2-4):\n")
print(sil_out)
cat("\nCluster composition at k = 3:\n")
print(counts)
cat("\nCross-region membership (head):\n")
print(head(cross, 10))
cat("\nFiles written to:", coa_dir, "\n")
