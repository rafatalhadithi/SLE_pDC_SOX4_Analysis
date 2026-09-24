# 09_generate_publication_figures.R
# Load required libraries
suppressPackageStartupMessages({
  library(reticulate)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(viridis)
})

# Force R to use the Python from your active environment
use_python(Sys.which("python"), required = TRUE)

# --- 1. SETUP AND DATA LOADING ---
cat("1. Initializing Python connection via reticulate...\n")
sc <- import("scanpy")
pd <- import("pandas")
ad <- import("anndata")

# Use dynamic relative paths (Run this from the repository root)
input_file <- "data/SLE_pDC_Integrated.h5ad"
output_dir <- "figures/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("2. Loading integrated dataset (this may take a moment)...\n")
adata <- sc$read_h5ad(input_file)

# --- 2. FIGURE 1A: UMAP VALIDATION ---
cat("3. Generating Figure 1A (UMAPs)...\n")
umap_coords <- adata$obsm[["X_pca_harmony_umap"]] # Or "X_umap" depending on scanpy version output
if (is.null(umap_coords)) {
  umap_coords <- adata$obsm[["X_umap"]]
}

umap_df <- data.frame(
  UMAP1 = umap_coords[,1],
  UMAP2 = umap_coords[,2]
)

features <- c("SOX4(+)", "CXCR4", "STAT1(+)", "ISG15")
for (feat in features) {
  umap_df[[feat]] <- as.numeric(adata$obs_vector(feat))
}

plot_umap <- function(df, feature_name, title) {
  df <- df %>% arrange(.data[[feature_name]])
  
  ggplot(df, aes(x = UMAP1, y = UMAP2, color = .data[[feature_name]])) +
    geom_point(size = 0.6, alpha = 0.85, stroke = 0) +
    scale_color_gradientn(
      colors = c("grey90", "#FDC926", "#FA9E3B", "#ED7953", "#D8576B", "#BD3786", "#9C179E", "#7000A8", "#46039F"),
      name = "Expr"
    ) + 
    theme_void(base_size = 14) + 
    labs(title = title) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15, margin = margin(b = 10)),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "right",
      legend.key.height = unit(1.2, "cm"),
      legend.key.width = unit(0.3, "cm"),
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      plot.margin = margin(10, 10, 10, 10)
    )
}

p1 <- plot_umap(umap_df, "SOX4(+)", "SOX4 Network Activity")
p2 <- plot_umap(umap_df, "CXCR4", "CXCR4 Gene Expression")
p3 <- plot_umap(umap_df, "STAT1(+)", "STAT1 Network Activity")
p4 <- plot_umap(umap_df, "ISG15", "ISG15 Gene Expression")

umap_grid <- (p1 | p2) / (p3 | p4) + plot_annotation(tag_levels = 'A')

ggsave(paste0(output_dir, "Fig_1A_Network_vs_Gene_Validation.pdf"), plot = umap_grid, width = 12, height = 10, dpi = 300)

# --- 3. FIGURE 1B: TOP REGULONS DOTPLOT ---
cat("4. Calculating Top Regulons for Figure 1B...\n")
py$adata <- adata

py_run_string("
import scanpy as sc
import pandas as pd
import anndata as ad

regulon_cols = [col for col in adata.obs.columns if col.endswith('(+)')]
auc_matrix = adata.obs[regulon_cols].copy()
clean_obs = adata.obs.drop(columns=regulon_cols)
adata_auc = ad.AnnData(X=auc_matrix, obs=clean_obs)

sc.tl.rank_genes_groups(
    adata_auc, 
    groupby='disease', 
    groups=['systemic lupus erythematosus'], 
    reference='normal', 
    method='wilcoxon'
)

top_regulons = pd.DataFrame(adata_auc.uns['rank_genes_groups']['names'])['systemic lupus erythematosus'].head(10).tolist()

records = []
for disease_status in ['normal', 'systemic lupus erythematosus']:
    subset_matrix = auc_matrix[clean_obs['disease'] == disease_status]
    for reg in top_regulons:
        expr_vals = subset_matrix[reg]
        records.append({
            'Regulon': reg,
            'Disease': 'Healthy' if disease_status == 'normal' else 'SLE',
            'Mean_Expr': expr_vals.mean(),
            'Pct_Expr': (expr_vals > 0).mean() * 100
        })
dot_df = pd.DataFrame(records)
")

dot_df <- py$dot_df
top_regulons <- py$top_regulons
dot_df$Regulon <- factor(dot_df$Regulon, levels = top_regulons)

cat("5. Plotting Figure 1B Dotplot...\n")
p_dot <- ggplot(dot_df, aes(x = Regulon, y = Disease)) +
  geom_point(aes(size = Pct_Expr, color = Mean_Expr)) +
  scale_color_gradientn(
    colors = c("grey95", "#FFE0B2", "#FF9800", "#E65100", "#BF360C"),
    name = "Mean Score"
  ) +
  scale_size_continuous(range = c(2.5, 9.5), name = "% Cells Active") + 
  theme_minimal(base_size = 14) +
  labs(x = "Top 10 SLE Regulon Networks", y = "") +
  theme(
    axis.text.x = element_text(face = "bold", size = 12, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold", size = 13, color = "black"),
    panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    plot.margin = margin(t = 10, r = 20, b = 10, l = 10)
  )

ggsave(paste0(output_dir, "Fig_1B_Top_SLE_Regulons_Dotplot_Landscape.pdf"), plot = p_dot, width = 10, height = 4.5, dpi = 300)
cat("SUCCESS: High-resolution PDFs have been generated in your figures folder!\n")
