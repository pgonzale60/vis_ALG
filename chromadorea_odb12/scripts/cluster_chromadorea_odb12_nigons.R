#!/usr/bin/env Rscript
# cluster_chromadorea_odb12_nigons.R
# Replicate Gonzalez de la Rosa et al. (2021) unsupervised clustering on chromadorea_odb12 BUSCOs

suppressPackageStartupMessages({
  renv_lib <- "/Users/pgonzale60/Documents/research/Auanema_rhodense/auanema_rhodense_chromosomal_genome_research/renv/library/macos/R-4.5/aarch64-apple-darwin20"
  if (dir.exists(renv_lib)) {
    .libPaths(c(renv_lib, .libPaths()))
  }
  library(cluster)
  library(tidyverse)
})

cat("=== 1. Loading 77 Assemblies Inventory & Allocations ===\n")
inv <- read.delim("data/synteny/tol_77_chromosomal_assemblies_inventory.tsv", sep = "\t", stringsAsFactors = FALSE)
cat("Loaded inventory with", nrow(inv), "assemblies.\n")

# Load chromadorea_odb12 allocations
alloc <- read.delim("data/synteny/chromadorea_odb12_77assemblies_allocations.tsv.gz", sep = "\t", stringsAsFactors = FALSE)
cat("Loaded", nrow(alloc), "chromadorea_odb12 chromosomal allocations.\n")

# Load legacy nematoda_odb10 nigon map
legacy_map <- read.delim("data/synteny/gene2Nigon_busco20200927.tsv.gz", sep = "\t", stringsAsFactors = FALSE)
colnames(legacy_map) <- c("legacy_og", "legacy_nigon")
cat("Loaded legacy Nigon map with", nrow(legacy_map), "genes.\n")

# Load nematoda_odb10 allocations for ground truth cross-mapping
alloc_n10 <- read.delim("data/synteny/nematoda_odb10_77assemblies_allocations.tsv.gz", sep = "\t", stringsAsFactors = FALSE)
cat("Loaded", nrow(alloc_n10), "nematoda_odb10 chromosomal allocations.\n")

# Filter balanced panel of 18 high-quality representative genera across Chromadorea
selected_genera <- c(
  "Caenorhabditis", "Oscheius", "Auanema", "Litoditis",
  "Koerneria",
  "Cooperia", "Dictyocaulus", "Trichostrongylus", "Ostertagia",
  "Acanthocheilonema", "Dirofilaria", "Thelazia", "Spirocerca",
  "Acrobeloides", "Panagrellus", "Panagrolaimus", "Turbatrix", "Sphaerularia"
)

# For each genus, select the assembly with highest single-copy BUSCO count
panel <- inv %>%
  mutate(genus = sub("_.*", "", Species)) %>%
  filter(genus %in% selected_genera) %>%
  group_by(genus) %>%
  slice_max(BUSCO_Complete_SC, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("\n=== 2. Selected Phylogenetically Balanced Panel (", nrow(panel), "Genera) ===\n")
for (i in seq_len(nrow(panel))) {
  cat(sprintf("  %-18s | %-30s | %-15s | %d chrs | SC=%d (%s)\n",
              panel$genus[i], panel$Species[i], panel$Accession[i],
              panel$n_chrs[i], as.integer(panel$BUSCO_Complete_SC[i]), panel$BUSCO_SC_Pct[i]))
}

# Filter allocations for panel species
panel_accs <- panel$Accession
alloc_panel <- alloc %>% filter(accession %in% panel_accs)

# Create locus-by-species matrix of chromosome allocations
cat("\n=== 3. Building Ortholog-by-Species Allocation Matrix ===\n")
locus_counts <- alloc_panel %>%
  group_by(busco_id) %>%
  summarise(n_taxa = n_distinct(accession), .groups = "drop")

# 80% presence threshold
n_panel <- nrow(panel)
min_taxa <- ceiling(0.80 * n_panel)
cat(sprintf("Total BUSCO loci in allocations: %d\n", nrow(locus_counts)))
cat(sprintf("Applying >= 80%% presence filter (>= %d of %d panel species)...\n", min_taxa, n_panel))

core_loci <- locus_counts %>% filter(n_taxa >= min_taxa) %>% pull(busco_id)
cat(sprintf("Core BUSCO loci passing 80%% filter: %d / 1126 (%.1f%%)\n",
            length(core_loci), length(core_loci) / 1126 * 100))

# Filter allocations to core loci
core_alloc <- alloc_panel %>% filter(busco_id %in% core_loci)

# Build binary indicator matrix
# Features are (accession : chromosome)
core_alloc <- core_alloc %>%
  mutate(feat = paste(accession, chromosome, sep = "::"))

feats <- sort(unique(core_alloc$feat))
mat <- matrix(0L, nrow = length(core_loci), ncol = length(feats),
              dimnames = list(core_loci, feats))

idx <- cbind(match(core_alloc$busco_id, core_loci), match(core_alloc$feat, feats))
mat[idx] <- 1L

cat(sprintf("Allocation matrix dimensions: %d loci x %d chromosome features.\n", nrow(mat), ncol(mat)))

# Compute pairwise Dice distance matrix:
# Dice similarity = 2 * (u . v) / (|u|_1 + |v|_1)
# Dice distance = 1 - Dice similarity
cat("\n=== 4. Computing Pairwise Dice Distance Matrix ===\n")
# Matrix multiplication to get intersection count A
A <- mat %*% t(mat)
# Row sums (|u|_1)
row_sums <- rowSums(mat)
# Denominator matrix: |u|_1 + |v|_1
denom <- outer(row_sums, row_sums, "+")
# Dice distance
dice_dist <- 1 - (2 * A) / denom
diag(dice_dist) <- 0
dice_dist_obj <- as.dist(dice_dist)

cat("Dice distance computed successfully.\n")

# Run PAM clustering across k = 2 .. 10 and compute average silhouette width
cat("\n=== 5. Evaluating Cluster Robustness (k = 2 .. 10) ===\n")
cat("k\tAvg_Silhouette_Width\n")
sil_scores <- numeric(10)
pam_models <- list()

for (k in 2:10) {
  p_fit <- pam(dice_dist_obj, k = k, diss = TRUE)
  sil_scores[k] <- p_fit$silinfo$avg.width
  pam_models[[k]] <- p_fit
  cat(sprintf("%d\t%.4f%s\n", k, sil_scores[k], ifelse(k == 7, "  <-- NIGON HYPOTHESIS (k=7)", "")))
}

best_k <- which.max(sil_scores[2:10]) + 1
cat(sprintf("\nOptimal cluster number by maximum average silhouette width: k = %d (silhouette = %.4f)\n",
            best_k, sil_scores[best_k]))

# Inspect k = 7 clustering
fit7 <- pam_models[[7]]
cluster7 <- fit7$clustering

cat("\n=== 6. Cluster Sizes for k = 7 ===\n")
print(table(cluster7))

# Save output object for inspection
saveRDS(list(
  panel = panel,
  core_loci = core_loci,
  dice_dist = dice_dist,
  sil_scores = sil_scores,
  fit7 = fit7
), "data/synteny/pam_k7_clustering_result.rds")

cat("\nResults saved to data/synteny/pam_k7_clustering_result.rds\n")
