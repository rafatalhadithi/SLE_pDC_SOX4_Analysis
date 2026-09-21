#!/usr/bin/env Rscript
# Script: 10_generate_publication_figures.R
# Purpose: High-resolution vector graphic composite showing integration patterns and expressions.

library(reticulate)
library(ggplot2)
library(dplyr)
library(patchwork)
library(viridis)

# Dynamic Paths
base_dir <- getwd()
data_dir <- file.path(base_dir, "data")
figures_dir <- file.path(base_dir, "figures")
input_file <- file.path(data_dir, "SLE_pDC_Integrated.h5ad")

if (!file.exists(input_file)) {
  stop("❌ ERROR: Integrated h5ad file not found in data directory.")
}

cat("1. Configuring Python environment via reticulate...\n")
python_path <- Sys.which("python")
if (nchar(python_path) > 0) {
  use_python(python_path, required = TRUE)
} else {
  try(use_condaenv("base", required = TRUE), silent = TRUE)
}

sc <- import("scanpy")
cat("2. Loading dataset...\n")
adata <- sc$read_h5ad(input_file)

plot_data <- data.frame(
  UMAP1 = adata$obsm[["X_umap"]][, 1],
  UMAP2 = adata$obsm[["X_umap"]][, 2],
  Donor = as.character(adata$obs$donor_id),
  Disease = as.character(adata$obs$disease),
  SOX4 = as.numeric(adata$obs_vector("SOX4"))
)

cat("3. Rendering High-Resolution Vectors...\n")
umap_theme <- theme_classic(base_size = 12) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14))

p1 <- ggplot(plot_data, aes(x=UMAP1, y=UMAP2, color=Donor)) +
  geom_point(data = plot_data[sample(nrow(plot_data)), ], alpha=0.7, size=0.6, stroke=0) +
  labs(title="Integration by Donor") + umap_theme + theme(legend.position="none")

disease_colors <- c("systemic lupus erythematosus" = "#8E44AD", "normal" = "#7F8C8D")
p2 <- ggplot(plot_data, aes(x=UMAP1, y=UMAP2, color=Disease)) +
  geom_point(data = plot_data[sample(nrow(plot_data)), ], alpha=0.7, size=0.6, stroke=0) +
  scale_color_manual(values=disease_colors) + labs(title="Disease State") + umap_theme

p3 <- ggplot(plot_data, aes(x=UMAP1, y=UMAP2, color=SOX4)) +
  geom_point(data = plot_data[order(plot_data$SOX4), ], alpha=0.8, size=0.6, stroke=0) +
  scale_color_viridis(option="viridis", name="Expr") + labs(title="SOX4 Expression") + umap_theme

final_plot <- p1 + p2 + p3 + plot_annotation(tag_levels='A') & theme(plot.tag=element_text(face='bold', size=16))

output_pdf <- file.path(figures_dir, "Supplementary_FigS1_Harmony_Integration.pdf")
ggsave(output_pdf, plot=final_plot, width=15, height=4.5, device=cairo_pdf)

cat("SUCCESS: All figures generated.\n")
