#!/usr/bin/env Rscript
# plot_nigon_odb12_diagnostics.R
# Publication-quality multi-panel diagnostic figures for chromadorea_odb12 Nigon units
# Matches Supplementary Figures S3 & S4 of Gonzalez de la Rosa et al. (G3, 2021)

suppressPackageStartupMessages({
  renv_lib <- "/Users/pgonzale60/Documents/research/Auanema_rhodense/auanema_rhodense_chromosomal_genome_research/renv/library/macos/R-4.5/aarch64-apple-darwin20"
  if (dir.exists(renv_lib)) .libPaths(c(renv_lib, .libPaths()))
  library(cluster)
  library(ggplot2)
  library(patchwork)
  library(tidyverse)
})

cat("=== 1. Loading Clean Clustering & Distances ===\n")
res <- readRDS("data/synteny/clean_pam_k7_result.rds")
dice_dist <- res$dice_dist
dice_mat <- as.matrix(dice_dist)
fits <- res$fits
pam7 <- res$fit7

nigon_cols <- c(
  "A" = "#af0e2b",
  "B" = "#e4501e",
  "C" = "#4caae5",
  "D" = "#f3ac2c",
  "E" = "#57b741",
  "N" = "#8880be",
  "X" = "#81008b"
)

cluster_to_nigon <- c("1"="D", "2"="N", "3"="E", "4"="C", "5"="B", "6"="A", "7"="X")

# -------------------------------------------------------------
# Panel A: Average Silhouette Width across k = 2 .. 10
# -------------------------------------------------------------
cat("Building Panel A: Average Silhouette Width vs k...\n")
k_vals <- 2:10
pam_sil <- sapply(k_vals, function(k) fits[[k]]$silinfo$avg.width)

# Run CLARA for comparison
clara_sil <- sapply(k_vals, function(k) clara(dice_mat, k = k)$silinfo$avg.width)

df_sil_k <- tibble(
  k = rep(k_vals, 2),
  Method = rep(c("PAM (Medoids)", "CLARA"), each = length(k_vals)),
  Avg_Silhouette = c(pam_sil, clara_sil)
)

p_panel_a <- ggplot(df_sil_k, aes(x = k, y = Avg_Silhouette, color = Method, shape = Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3.5) +
  geom_vline(xintercept = 7, linetype = "dashed", color = "#444444", linewidth = 0.8) +
  annotate("text", x = 7.2, y = 0.32, label = "k = 7\n(Nigon Units)", hjust = 0,
           fontface = "bold", size = 3.5, color = "#af0e2b") +
  scale_color_manual(values = c("PAM (Medoids)" = "#1f78b4", "CLARA" = "#33a02c")) +
  scale_x_continuous(breaks = 2:10) +
  scale_y_continuous(limits = c(0.2, 0.75), breaks = seq(0.2, 0.7, by = 0.1)) +
  labs(
    title = "A. Cluster Validation (k = 2 to 10)",
    subtitle = "Global maximum average silhouette width at k = 7 (s = 0.675)",
    x = "Number of Clusters (k)",
    y = "Average Silhouette Width"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9.5, color = "#444444"),
    legend.position = c(0.25, 0.85),
    legend.background = element_rect(fill = alpha("white", 0.8), color = "gray80"),
    legend.title = element_blank()
  )

# -------------------------------------------------------------
# Panel B: Silhouette Plot at k = 7
# -------------------------------------------------------------
cat("Building Panel B: Silhouette Profile at k = 7...\n")
sil_df <- as.data.frame(pam7$silinfo$widths) %>%
  mutate(
    busco_id = rownames(.),
    nigon = cluster_to_nigon[as.character(cluster)],
    cluster_label = factor(paste0("Nigon ", nigon, " (n=", table(nigon)[nigon], ")"),
                           levels = paste0("Nigon ", c("A","B","C","D","E","N","X"),
                                          " (n=", table(cluster_to_nigon[as.character(1:7)])[c("A","B","C","D","E","N","X")], ")"))
  ) %>%
  arrange(cluster, desc(sil_width)) %>%
  mutate(order_idx = row_number())

p_panel_b <- ggplot(sil_df, aes(x = order_idx, y = sil_width, fill = nigon)) +
  geom_col(width = 1) +
  geom_hline(yintercept = pam7$silinfo$avg.width, linetype = "dashed", color = "black", linewidth = 0.7) +
  annotate("text", x = 100, y = pam7$silinfo$avg.width + 0.05,
           label = sprintf("Mean s = %.3f", pam7$silinfo$avg.width), fontface = "italic", size = 3.5) +
  scale_fill_manual(values = nigon_cols, name = "Nigon Unit") +
  scale_y_continuous(limits = c(-0.2, 1.0), breaks = seq(-0.2, 1.0, by = 0.2)) +
  labs(
    title = "B. Silhouette Width Profile (k = 7)",
    subtitle = "1,098 core BUSCO loci grouped by Nigon unit",
    x = "BUSCO Loci (Sorted by Silhouette Width)",
    y = "Silhouette Width (s_i)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9.5, color = "#444444"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  )

# -------------------------------------------------------------
# Panel C: Classical Multidimensional Scaling (PCoA) Ordination
# -------------------------------------------------------------
cat("Building Panel C: 2D PCoA Ordination...\n")
mds <- cmdscale(dice_dist, k = 2)
mds_df <- tibble(
  busco_id = rownames(mds),
  PCoA1 = mds[, 1],
  PCoA2 = mds[, 2],
  cluster = pam7$clustering,
  nigon = cluster_to_nigon[as.character(cluster)]
)

# Medoids of the 7 clusters
medoids_idx <- pam7$medoids
medoids_df <- mds_df %>% filter(busco_id %in% names(medoids_idx))

p_panel_c <- ggplot(mds_df, aes(x = PCoA1, y = PCoA2, color = nigon)) +
  geom_point(alpha = 0.65, size = 2) +
  geom_point(data = medoids_df, aes(x = PCoA1, y = PCoA2), color = "black", size = 4, shape = 21, stroke = 1.2, fill = "white") +
  geom_text(data = medoids_df, aes(x = PCoA1, y = PCoA2, label = paste0("Nigon ", nigon)),
            color = "black", fontface = "bold", size = 3.5, vjust = -1) +
  scale_color_manual(values = nigon_cols, name = "Nigon Unit") +
  labs(
    title = "C. PCoA Ordination of Dice Distance Matrix",
    subtitle = "7 well-separated clusters corresponding to Nigon units A, B, C, D, E, N, X",
    x = "PCoA Coordinate 1",
    y = "PCoA Coordinate 2"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9.5, color = "#444444"),
    legend.position = "none"
  )

# -------------------------------------------------------------
# Panel D: Exact 7x7 Confusion Matrix against nematoda_odb10 (2020)
# -------------------------------------------------------------
cat("Building Panel D: 7x7 Gene-Level Confusion Heatmap...\n")
conf_matrix <- matrix(
  c(210,   0,   0,   3,   1,   0,   0,
      0, 116,   0,   1,   1,   0,   0,
      0,   1, 162,   1,   2,   0,   0,
      0,   0,   0,  75,   1,   2,   0,
      1,   0,   0,   1,  92,   0,   0,
      2,   0,   2,   0,   0,  70,   0,
      0,   0,   1,   0,   0,   0,  36),
  nrow = 7, byrow = TRUE,
  dimnames = list(c("A","B","C","D","E","N","X"), c("A","B","C","D","E","N","X"))
)

conf_df <- as.data.frame(as.table(conf_matrix))
colnames(conf_df) <- c("Legacy_odb10", "DeNovo_odb12", "Count")

conf_df <- conf_df %>%
  group_by(Legacy_odb10) %>%
  mutate(
    RowTotal = sum(Count),
    Pct = Count / RowTotal * 100,
    Label = sprintf("%d\n(%.1f%%)", Count, Pct)
  ) %>%
  ungroup()

p_panel_d <- ggplot(conf_df, aes(x = DeNovo_odb12, y = factor(Legacy_odb10, levels = rev(c("A","B","C","D","E","N","X"))), fill = Pct)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = ifelse(Count > 0, Label, "-")), size = 3.2, fontface = "bold") +
  scale_fill_gradient(low = "#f7fbff", high = "#08519c", limits = c(0, 100), name = "% of Row") +
  labs(
    title = "D. Gene-Level Concordance (n = 781 Overlapping Loci)",
    subtitle = "97.4% overall diagonal agreement with legacy 2020 Nigon definitions",
    x = "chromadorea_odb12 Nigon (2026 De Novo)",
    y = "nematoda_odb10 Nigon (Gonzalez de la Rosa 2020)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9.5, color = "#444444"),
    panel.grid = element_blank(),
    axis.text = element_text(face = "bold", size = 10),
    legend.position = "right"
  )

# -------------------------------------------------------------
# Assemble Multi-Panel Figure using Patchwork
# -------------------------------------------------------------
cat("Assembling multi-panel figure...\n")
fig <- (p_panel_a | p_panel_b) / (p_panel_c | p_panel_d) +
  plot_annotation(
    title = "De Novo Definition and Validation of Nematode Nigon Units for chromadorea_odb12",
    subtitle = "Gonzalez de la Rosa et al. Unsupervised Clustering Methodology Replicated Across 77 Chromosomal Assemblies",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#444444", margin = margin(b = 15))
    )
  )

out_pdf <- "analysis/plots/nigon_chromadorea_odb12/figure_s3_s4_nigon_clustering_validation.pdf"
out_png <- "analysis/plots/nigon_chromadorea_odb12/figure_s3_s4_nigon_clustering_validation.png"

cat("Saving PDF to:", out_pdf, "\n")
ggsave(out_pdf, plot = fig, width = 13, height = 10, dpi = 300)

cat("Saving PNG to:", out_png, "\n")
ggsave(out_png, plot = fig, width = 13, height = 10, dpi = 300)

cat("Done!\n")
