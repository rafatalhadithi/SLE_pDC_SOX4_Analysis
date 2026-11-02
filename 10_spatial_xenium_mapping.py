import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu
import os
import warnings

# Suppress annoying seaborn deprecation warnings
warnings.filterwarnings("ignore")

print("--- Step 1: Initializing Spatial Analysis Environment ---")
data_path = "F:/Lupus/GSE294965_Spatial/GSE294965_processed_data.h5ad"
output_dir = "E:/OneDrive/Review/Research idea/Bioinformatics/Lupus/reviewer_revisions/figures"
os.makedirs(output_dir, exist_ok=True)

print(f"Loading high-resolution 10x Xenium spatial matrix from:\n{data_path}...")
adata = sc.read_h5ad(data_path)

print("\n--- Step 2: Extracting Target Metadata ---")
cell_type_col = 'celltype_l1'
disease_col = 'Disease'
sample_col = 'sample' if 'sample' in adata.obs.columns else 'Slide_ID'

ln_mask = adata.obs[disease_col].astype(str).str.contains('LN|SLE|Lupus', case=False, na=False)
adata_ln = adata[ln_mask].copy() if ln_mask.sum() > 0 else adata.copy()

print(f"Total Lupus Nephritis cells across all biopsies: {adata_ln.n_obs}")
print(f"Author's exact cell type labels: {adata_ln.obs[cell_type_col].unique().tolist()}")

print("\n--- Step 3: Aggregating Spatial Distances Across ALL Biopsies ---")
pdc_distances = []
bg_distances = []

unique_samples = adata_ln.obs[sample_col].unique()
print(f"Looping through {len(unique_samples)} LN tissue cores for statistical power...")

best_sample = unique_samples[0]
max_combined_cells = 0

for sample in unique_samples:
    adata_sub = adata_ln[adata_ln.obs[sample_col] == sample]
    
    # 1. Identify pDCs
    pdc_mask = adata_sub.obs[cell_type_col].astype(str).str.contains('pDC|Plasmacytoid', case=False, na=False)
    
    # 2. Identify Target Niches using the proven CXCL12 gene expression
    if 'CXCL12' in adata_sub.var_names:
        cxcl12_expr = adata_sub[:, 'CXCL12'].X
        if hasattr(cxcl12_expr, 'toarray'): cxcl12_expr = cxcl12_expr.toarray().flatten()
        else: cxcl12_expr = cxcl12_expr.flatten()
        niche_mask = cxcl12_expr > 0
    else:
        # Fallback if CXCL12 is somehow missing (it shouldn't be based on previous run)
        niche_mask = adata_sub.obs[cell_type_col].astype(str).str.contains('Endo|Strom|Vessel|Fibro', case=False, na=False)
    
    bg_mask = ~(pdc_mask | niche_mask)
    
    # Track the best sample for plotting Panel A
    combined_cells = pdc_mask.sum() + niche_mask.sum()
    if combined_cells > max_combined_cells and pdc_mask.sum() > 0 and niche_mask.sum() > 0:
        max_combined_cells = combined_cells
        best_sample = sample
    
    # Skip cores that lack either pDCs or CXCL12 cells to prevent math errors
    if pdc_mask.sum() == 0 or niche_mask.sum() == 0:
        continue
        
    coords = adata_sub.obs[['x_centroid', 'y_centroid']].values
    pdc_coords = coords[pdc_mask]
    niche_coords = coords[niche_mask]
    bg_coords = coords[bg_mask]
    
    # Calculate nearest-neighbor distances
    tree = cKDTree(niche_coords)
    d_pdc, _ = tree.query(pdc_coords)
    pdc_distances.extend(d_pdc)
    
    d_bg, _ = tree.query(bg_coords)
    bg_distances.extend(d_bg)

print(f"\nAggregated N = {len(pdc_distances)} pDCs and N = {len(bg_distances)} Background Cells.")
stat, p_val = mannwhitneyu(pdc_distances, bg_distances, alternative='less')
print(f"Mann-Whitney U Test p-value: {p_val}")

print(f"\n--- Step 4: Generating 2-Panel Publication Figure (Using Best Sample: {best_sample}) ---")
plt.rcParams.update({'font.sans-serif': 'Arial'})
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7), dpi=300, gridspec_kw={'width_ratios': [1.2, 1]})

# PANEL A: Representative Spatial Map
adata_plot = adata_ln[adata_ln.obs[sample_col] == best_sample]

pdc_plot_mask = adata_plot.obs[cell_type_col].astype(str).str.contains('pDC|Plasmacytoid', case=False, na=False)
if 'CXCL12' in adata_plot.var_names:
    cxcl12_expr = adata_plot[:, 'CXCL12'].X
    if hasattr(cxcl12_expr, 'toarray'): cxcl12_expr = cxcl12_expr.toarray().flatten()
    else: cxcl12_expr = cxcl12_expr.flatten()
    niche_plot_mask = cxcl12_expr > 0
else:
    niche_plot_mask = adata_plot.obs[cell_type_col].astype(str).str.contains('Endo|Strom|Vessel|Fibro', case=False, na=False)

bg_plot_mask = ~(pdc_plot_mask | niche_plot_mask)

x_c = adata_plot.obs['x_centroid'].values
y_c = adata_plot.obs['y_centroid'].values

ax1.scatter(x_c[bg_plot_mask], y_c[bg_plot_mask], c='#e0e0e0', s=15, alpha=0.8, label='Kidney Tissue')
ax1.scatter(x_c[niche_plot_mask], y_c[niche_plot_mask], c='#1f77b4', s=40, alpha=0.9, label='CXCL12-Secreting Niche')
ax1.scatter(x_c[pdc_plot_mask], y_c[pdc_plot_mask], c='#d62728', s=120, alpha=1.0, label='Infiltrating pDCs', edgecolors='white', linewidth=1.2)

ax1.set_aspect('equal')
ax1.axis('off')
ax1.set_title(f"A. Representative pDC Infiltration in Lupus Nephritis\n(10x Xenium | TMA Core: {best_sample})", fontsize=14, fontweight='bold', pad=10)
ax1.legend(loc='lower center', bbox_to_anchor=(0.5, -0.15), scatterpoints=1, frameon=False, markerscale=1.5, fontsize=12, ncol=3)

# PANEL B: Statistical Aggregation Boxplot
df_stats = pd.DataFrame({
    'Distance (\u03bcm)': np.concatenate([pdc_distances, bg_distances]),
    'Cell Type': ['Infiltrating pDCs'] * len(pdc_distances) + ['Random Kidney Tissue'] * len(bg_distances)
})

sns.boxplot(data=df_stats, x='Cell Type', y='Distance (\u03bcm)', ax=ax2, hue='Cell Type', palette=['#d62728', '#e0e0e0'], legend=False, showfliers=False)
ax2.set_title("B. Pan-Cohort Spatial Colocalization\n(Aggregated Nearest-Neighbor Distances)", fontsize=14, fontweight='bold', pad=10)
ax2.set_ylabel("Distance to Nearest CXCL12 Niche", fontsize=12, fontweight='bold')
ax2.set_xlabel("")
ax2.tick_params(labelsize=12)

def format_pval(p):
    return f"p = {p:.1e}" if p < 0.001 else f"p = {p:.3f}"
    
ax2.text(0.5, ax2.get_ylim()[1] * 0.95, format_pval(p_val), ha='center', va='top', fontsize=14, fontweight='bold')

plt.tight_layout()

output_pdf = os.path.normpath(f"{output_dir}/Figure_6_Spatial_Transcriptomics_Full.pdf")
output_png = os.path.normpath(f"{output_dir}/Figure_6_Spatial_Transcriptomics_Full.png")

plt.savefig(output_pdf, format='pdf', bbox_inches='tight')
plt.savefig(output_png, dpi=300, bbox_inches='tight')
plt.close()

print(f"\nSUCCESS: 2-Panel Q1 Figure generated!")
print(f"Saved to: {output_pdf}")