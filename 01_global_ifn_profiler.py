#!/usr/bin/env python3
"""
Script: 01_global_ifn_profiler.py
Purpose: Evaluates the systemic contribution of various immune subsets to the circulating 
         Type I IFN pool using the 1.2M cell global SLE atlas (Perez et al. 2022).
"""

import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys

# -------------------------------------------------------------------------
# 1. Dynamic Path Configuration (Reviewer-Safe)
# -------------------------------------------------------------------------
# Automatically build paths relative to where the script is executed
BASE_DIR = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
DATA_DIR = os.path.join(BASE_DIR, "data")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")
FIGURES_DIR = os.path.join(BASE_DIR, "figures")

# Auto-create output directories
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

global_h5ad_path = os.path.join(DATA_DIR, "4118e166-34f5-4c1f-9eed-c64b90a3dace.h5ad")

if not os.path.exists(global_h5ad_path):
    print(f"❌ ERROR: Global atlas not found at {global_h5ad_path}")
    print("Please download the CZ CELLxGENE dataset and place it in the 'data' folder.")
    sys.exit(1)

# -------------------------------------------------------------------------
# 2. Reading Global Atlas Metadata & Scanning IFN Genes
# -------------------------------------------------------------------------
print("--- Step 1: Reading Global Atlas Metadata & Scanning IFN Genes ---")
adata_backed = sc.read_h5ad(global_h5ad_path, backed='r')

canonical_ifn_genes = [
    'IFNA1', 'IFNA2', 'IFNA4', 'IFNA5', 'IFNA6', 'IFNA7', 
    'IFNA8', 'IFNA10', 'IFNA13', 'IFNA14', 'IFNA16', 'IFNA17', 
    'IFNA21', 'IFNB1', 'IFNW1'
]

# Identify where gene symbols are stored and map to Ensembl IDs if necessary
if 'feature_name' in adata_backed.var.columns:
    mask = adata_backed.var['feature_name'].isin(canonical_ifn_genes)
    var_names_to_fetch = adata_backed.var.index[mask].tolist()
elif 'gene_symbol' in adata_backed.var.columns:
    mask = adata_backed.var['gene_symbol'].isin(canonical_ifn_genes)
    var_names_to_fetch = adata_backed.var.index[mask].tolist()
else:
    var_names_to_fetch = [g for g in canonical_ifn_genes if g in adata_backed.var_names]

print(f"Found {len(var_names_to_fetch)} Type I IFN genes mapped to index IDs.")

# -------------------------------------------------------------------------
# 3. Extracting IFN Expression Matrix
# -------------------------------------------------------------------------
print("--- Step 2: Extracting Metadata and IFN Expression Matrix ---")
cell_type_col = 'cell_type' if 'cell_type' in adata_backed.obs.columns else 'author_cell_type'
disease_col = 'disease' if 'disease' in adata_backed.obs.columns else 'disease_state'

df_meta = adata_backed.obs[[cell_type_col, disease_col]].copy()
df_meta.columns = ['Cell_Type', 'Disease']

# Bypass the AnnData backed .raw bug by only extracting the numeric matrix
if adata_backed.raw is not None:
    print("Extracting numeric counts directly from adata.raw.X...")
    expr_matrix = adata_backed.raw[:, var_names_to_fetch].X
else:
    print("adata.raw not found. Extracting from adata.X...")
    expr_matrix = adata_backed[:, var_names_to_fetch].X

# Safely calculate total Type I IFN expression per cell
if sp.issparse(expr_matrix):
    raw_sums = np.asarray(expr_matrix.sum(axis=1)).squeeze()
else:
    raw_sums = np.asarray(expr_matrix.sum(axis=1)).squeeze()

df_meta['Total_IFN_Expr'] = raw_sums
print(f"Total Type I IFN transcripts detected: {raw_sums.sum():.2f}")

# -------------------------------------------------------------------------
# 4. Statistical Breakdown by Cell Type in SLE
# -------------------------------------------------------------------------
print("\n--- Step 3: Statistical Breakdown by Cell Type in SLE ---")
sle_mask = df_meta['Disease'].astype(str).str.contains('lupus|SLE|Systemic Lupus', case=False, na=False)
if sle_mask.sum() > 0:
    df_sle = df_meta[sle_mask].copy()
    print(f"Analyzed {len(df_sle):,} cells from SLE patients.")
else:
    df_sle = df_meta.copy()

# Group by Cell Type
grouped = df_sle.groupby('Cell_Type', observed=False)['Total_IFN_Expr'].agg(['sum', 'mean', 'count'])
grouped = grouped[grouped['count'] > 0] 

total_ifn_sle = grouped['sum'].sum()
if total_ifn_sle > 0:
    grouped['Percent_Share'] = (grouped['sum'] / total_ifn_sle) * 100
else:
    grouped['Percent_Share'] = 0.0

grouped = grouped.sort_values(by='Percent_Share', ascending=False)

print("\n=======================================================")
print(" TYPE I IFN PRODUCTION SHARE BY CELL TYPE (SLE)      ")
print("=======================================================")
print(grouped[['count', 'mean', 'Percent_Share']].to_string())
print("=======================================================\n")

csv_out = os.path.join(OUTPUT_DIR, "01_Global_Type1_IFN_Production_Summary.csv")
grouped.to_csv(csv_out)

# -------------------------------------------------------------------------
# 5. Generating Publication-Quality Figure
# -------------------------------------------------------------------------
print("--- Step 4: Generating Publication-Quality Figure ---")
if total_ifn_sle > 0:
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    plt.rcParams.update({'font.size': 11, 'font.sans-serif': 'Arial'})

    top10_share = grouped.head(10).reset_index()

    sns.barplot(data=top10_share, x='Percent_Share', y='Cell_Type', hue='Cell_Type', dodge=False, width=0.7, ax=axes[0], palette='Blues_r', legend=False)
    axes[0].set_title('A. Share of Total Type I IFN Transcripts in SLE (%)', fontweight='bold', fontsize=12)
    axes[0].set_xlabel('Percentage Contribution (%)')
    axes[0].set_ylabel('Cell Type')

    for i in range(len(top10_share)):
        val = top10_share.iloc[i]['Percent_Share']
        if pd.notna(val):
            axes[0].text(val + 0.5, i, f"{val:.1f}%", va='center', fontsize=10)

    sns.barplot(data=top10_share, x='mean', y='Cell_Type', hue='Cell_Type', dodge=False, width=0.7, ax=axes[1], palette='Reds_r', legend=False)
    axes[1].set_title('B. Mean Type I IFN Expression per Cell', fontweight='bold', fontsize=12)
    axes[1].set_xlabel('Mean Normalized Expression')
    axes[1].set_ylabel('')

    plt.tight_layout()
    
    fig_out_png = os.path.join(FIGURES_DIR, "01_Global_Type1_IFN_Production.png")
    fig_out_pdf = os.path.join(FIGURES_DIR, "01_Global_Type1_IFN_Production.pdf")
    
    plt.savefig(fig_out_png, dpi=300, bbox_inches='tight')
    plt.savefig(fig_out_pdf, format='pdf', bbox_inches='tight')
    plt.close()
    
    print(f"Success! Visuals saved to: {FIGURES_DIR}")
else:
    print("Warning: No IFN expression detected, skipping plot generation.")
