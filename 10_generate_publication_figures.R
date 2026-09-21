#!/usr/bin/env Rscript
# Script: 10_generate_publication_figures.R
# Purpose: High-resolution vector graphic composite showing integration patterns and expressions.

suppressPackageStartupMessages({
  library(reticulate)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(viridis)
})

# Dynamic Paths
base_dir <- getwd()
data_dir <- file.path(base_dir, "data")
figures_dir <- file.path(base_dir, "figures")
dir.create(figures_dir, showWarnings = FALSE)

input_integrated <- file.path(data_dir, "SLE_pDC_Integrated.h5ad")
input_scenic <- file.path(data_dir, "SLE_pDC_pySCENIC_Complete.h5ad")

if (!file.exists(input_integrated) | !file.exists(input_scenic)) {
  stop("❌ ERROR: Required .h5ad files not found in data/ directory.")
}

cat("1. Configuring Python environment via reticulate...\n")
python_path <- Sys.which("python")
if (nchar(python_path) > 0) {
  use_python(python_path, required = TRUE)
} else {
  try(use_condaenv("base", required = TRUE), silent = TRUE)
}

sc <- import("scanpy")
ad <- import("anndata")
pd <- import("pandas")

# ==============================================================================
# PART 1: SUPPLEMENTARY FIGURE S1 (Harmony Integration)
# ==============================================================================
cat("2. Generating Supplementary Fig S1 (Harmony UMAP)...\n")
adata_int <- sc$read_h5ad(input_integrated)

plot_data_int <- data.frame(
  UMAP1 = adata_int$obsm[["X_umap"]][, 1],
  UMAP2 = adata_int$obsm[["X_umap"]][, 2],
  Donor = as.character(adata_int$obs$donor_id),
  Disease = as.character(adata_int$obs$disease),
  SOX4 = as.numeric(adata_int$obs_vector("SOX4"))
)

umap_theme <- theme_classic(base_size = 12) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14))

p1 <- ggplot(plot_data_int, aes(x=UMAP1, y=UMAP2, color=Donor)) +
  geom_point(data = plot_data_int[sample(nrow(plot_data_int)), ], alpha=0.7, size=0.6, stroke=0) +
  labs(title="Integration by Donor") + umap_theme + theme(legend.position="none")

disease_colors <- c("systemic lupus erythematosus" = "#8E44AD", "normal" = "#7F8C8D")
p2 <- ggplot(plot_data_int, aes(x=UMAP1, y=UMAP2, color=Disease)) +
  geom_point(data = plot_data_int[sample(nrow(plot_data_int)), ], alpha=0.7, size=0.6, stroke=0) +
  scale_color_manual(values=disease_colors) + labs(title="Disease State") + umap_theme

p3 <- ggplot(plot_data_int, aes(x=UMAP1, y=UMAP2, color=SOX4)) +
  geom_point(data = plot_data_int[order(plot_data_int$SOX4), ], alpha=0.8, size=0.6, stroke=0) +
  scale_color_viridis(option="viridis", name="Expr") + labs(title="SOX4 Expression") + umap_theme

fig_s1 <- p1 + p2 + p3 + plot_annotation(tag_levels='A') & theme(plot.tag=element_text(face='bold', size=16))
ggsave(file.path(figures_dir, "Supplementary_FigS1_Harmony_Integration.pdf"), plot=fig_s1, width=15, height=4.5, device=cairo_pdf)

# ==============================================================================
# PART 2: FIGURE 1A & 1B (Network Validation & Dotplot)
# ==============================================================================
cat("3. Generating Figure 1 (Network Validation & Dotplot)...\n")
adata_scenic <- sc$read_h5ad(input_scenic)

# Extract UMAP coordinates safely
umap_df <- data.frame(
  UMAP1 = adata_scenic$obsm[["X_umap"]][, 1],
  UMAP2 = adata_scenic$obsm[["X_umap"]][, 2],
  SOX4_reg = as.numeric(adata_scenic$obs_vector("SOX4(+)")),
  CXCR4 = as.numeric(adata_scenic$obs_vector("CXCR4")),
  STAT1_reg = as.numeric(adata_scenic$obs_vector("STAT1(+)")),
  ISG15 = as.numeric(adata_scenic$obs_vector("ISG15"))
)

plot_umap_feat <- function(df, feature_name, title) {
  df <- df %>% arrange(.data[[feature_name]])
  ggplot(df, aes(x = UMAP1, y = UMAP2, color = .data[[feature_name]])) +
    geom_point(size = 0.6, alpha = 0.85, stroke = 0) +
    scale_color_gradientn(colors = c("grey90", "#FDC926", "#FA9E3B", "#ED7953", "#D8576B", "#BD3786", "#9C179E", "#7000A8", "#46039F"), name = "Expr") + 
    theme_void(base_size = 14) + labs(title = title) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15, margin = margin(b = 10)),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.position = "right", legend.key.height = unit(1.2, "cm"))
}

p_sox4 <- plot_umap_feat(umap_df, "SOX4_reg", "SOX4(+) Network Activity")
p_cxcr4 <- plot_umap_feat(umap_df, "CXCR4", "CXCR4 Gene Expression")
p_stat1 <- plot_umap_feat(umap_df, "STAT1_reg", "STAT1(+) Network Activity")
p_isg15 <- plot_umap_feat(umap_df, "ISG15", "ISG15 Gene Expression")

fig_1a <- (p_sox4 | p_cxcr4) / (p_stat1 | p_isg15) + plot_annotation(tag_levels = 'A')
ggsave(file.path(figures_dir, "Fig1A_Network_vs_Gene_Validation.pdf"), plot = fig_1a, width = 12, height = 10, device=cairo_pdf)

# --- DOTPLOT NATIVE PYTHON MATH ---
py$adata_scenic <- adata_scenic
py_run_string("
import scanpy as sc
import pandas as pd
import anndata as ad
regulon_cols = [col for col in adata_scenic.obs.columns if col.endswith('(+)')]
auc_matrix = adata_scenic.obs[regulon_cols].copy()
clean_obs = adata_scenic.obs.drop(columns=regulon_cols)
adata_auc = ad.AnnData(X=auc_matrix, obs=clean_obs)
sc.tl.rank_genes_groups(adata_auc, groupby='disease', groups=['systemic lupus erythematosus'], reference='normal', method='wilcoxon')
top_regulons = pd.DataFrame(adata_auc.uns['rank_genes_groups']['names'])['systemic lupus erythematosus'].head(10).tolist()
records = []
for disease_status in ['normal', 'systemic lupus erythematosus']:
    subset_matrix = auc_matrix[clean_obs['disease'] == disease_status]
    for reg in top_regulons:
        expr_vals = subset_matrix[reg]
        records.append({'Regulon': reg, 'Disease': 'Healthy' if disease_status == 'normal' else 'SLE', 'Mean_Expr': expr_vals.mean(), 'Pct_Expr': (expr_vals > 0).mean() * 100})
dot_df = pd.DataFrame(records)
")

dot_df <- py$dot_df
top_regulons <- py$top_regulons
dot_df$Regulon <- factor(dot_df$Regulon, levels = top_regulons)

fig_1b <- ggplot(dot_df, aes(x = Regulon, y = Disease)) +
  geom_point(aes(size = Pct_Expr, color = Mean_Expr)) +
  scale_color_gradientn(colors = c("grey95", "#FFE0B2", "#FF9800", "#E65100", "#BF360C"), name = "Mean Score") +
  scale_size_continuous(range = c(2.5, 9.5), name = "% Cells Active") + 
  theme_minimal(base_size = 14) + labs(x = "Top 10 SLE Regulon Networks", y = "") +
  theme(axis.text.x = element_text(face = "bold", size = 12, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(face = "bold", size = 13, color = "black"),
        panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "bottom")

ggsave(file.path(figures_dir, "Fig1B_Top_SLE_Regulons_Dotplot.pdf"), plot = fig_1b, width = 10, height = 4.5, device=cairo_pdf)

cat("SUCCESS: All publication vectors generated successfully.\n")
