#!/usr/bin/env python3
import scanpy as sc, pandas as pd, numpy as np, scipy.sparse as sp, matplotlib.pyplot as plt, seaborn as sns, os, sys

BASE_DIR = os.getcwd()
DATA_DIR, OUTPUT_DIR, FIGURES_DIR = [os.path.join(BASE_DIR, d) for d in ["data", "outputs", "figures"]]
for d in [DATA_DIR, OUTPUT_DIR, FIGURES_DIR]: os.makedirs(d, exist_ok=True)

global_h5ad_path = os.path.join(DATA_DIR, "4118e166-34f5-4c1f-9eed-c64b90a3dace.h5ad")

if __name__ == "__main__":
    print("1. Reading Global Atlas Metadata...")
    adata_backed = sc.read_h5ad(global_h5ad_path, backed='r')
    
    # --- PART A: GLOBAL IFN PROFILING ---
    ifn_genes = ['IFNA1', 'IFNA2', 'IFNA4', 'IFNA5', 'IFNA6', 'IFNA7', 'IFNA8', 'IFNA10', 'IFNA13', 'IFNA14', 'IFNA16', 'IFNA17', 'IFNA21', 'IFNB1', 'IFNW1']
    var_cols = adata_backed.var.columns
    name_col = 'feature_name' if 'feature_name' in var_cols else 'gene_symbol' if 'gene_symbol' in var_cols else None
    
    if name_col:
        mask = adata_backed.var[name_col].isin(ifn_genes)
        fetch_ids = adata_backed.var.index[mask].tolist()
    else:
        fetch_ids = [g for g in ifn_genes if g in adata_backed.var_names]

    df_meta = adata_backed.obs[['cell_type' if 'cell_type' in adata_backed.obs.columns else 'author_cell_type', 
                                'disease' if 'disease' in adata_backed.obs.columns else 'disease_state']].copy()
    df_meta.columns = ['Cell_Type', 'Disease']
    
    expr_matrix = adata_backed.raw[:, fetch_ids].X if adata_backed.raw is not None else adata_backed[:, fetch_ids].X
    df_meta['Total_IFN_Expr'] = np.asarray(expr_matrix.sum(axis=1)).squeeze() if sp.issparse(expr_matrix) else np.asarray(expr_matrix.sum(axis=1)).squeeze()

    sle_mask = df_meta['Disease'].astype(str).str.contains('lupus|SLE', case=False, na=False)
    grouped = df_meta[sle_mask].groupby('Cell_Type', observed=False)['Total_IFN_Expr'].agg(['sum', 'mean', 'count'])
    grouped = grouped[grouped['count'] > 0]
    grouped['Percent_Share'] = (grouped['sum'] / grouped['sum'].sum()) * 100 if grouped['sum'].sum() > 0 else 0.0
    grouped = grouped.sort_values(by='Percent_Share', ascending=False)
    
    grouped.to_csv(os.path.join(OUTPUT_DIR, "01_Global_Type1_IFN_Production_Summary.csv"))
    
    if grouped['sum'].sum() > 0:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        top10 = grouped.head(10).reset_index()
        sns.barplot(data=top10, x='Percent_Share', y='Cell_Type', hue='Cell_Type', ax=axes[0], palette='Blues_r', legend=False)
        axes[0].set_title('A. Share of Total Type I IFN Transcripts in SLE (%)', fontweight='bold')
        sns.barplot(data=top10, x='mean', y='Cell_Type', hue='Cell_Type', ax=axes[1], palette='Reds_r', legend=False)
        axes[1].set_title('B. Mean Type I IFN Expression per Cell', fontweight='bold')
        plt.tight_layout()
        plt.savefig(os.path.join(FIGURES_DIR, "01_Global_Type1_IFN_Production.pdf"), format='pdf', bbox_inches='tight')

    # --- PART B: PDC ISOLATION & PREP ---
    print("2. Isolating pDCs & Prepping for pySCENIC...")
    adata_pdc = adata_backed[(adata_backed.obs['sex'] == 'female') & (adata_backed.obs['cell_type'].str.contains('plasmacytoid', case=False, na=False))].to_memory()
    
    if name_col:
        adata_pdc.var_names = adata_pdc.var[name_col].values
    adata_pdc.var_names_make_unique()
    if adata_pdc.raw is not None: del adata_pdc.raw
    
    sc.pp.filter_cells(adata_pdc, min_counts=1)
    sc.pp.normalize_total(adata_pdc, target_sum=1e4)
    sc.pp.log1p(adata_pdc)
    
    sc.tl.rank_genes_groups(adata_pdc, groupby='disease', groups=['systemic lupus erythematosus'], reference='normal', method='wilcoxon')
    top_genes = sc.get.rank_genes_groups_df(adata_pdc, group='systemic lupus erythematosus')
    top77 = top_genes[(top_genes['logfoldchanges'] > 0.5) & (top_genes['pvals_adj'] < 0.05)]
    top77.to_csv(os.path.join(OUTPUT_DIR, "SLE_pDC_Top77_Upregulated.csv"), index=False)
    
    adata_pdc.var.index.name = None
    adata_pdc.write_h5ad(os.path.join(DATA_DIR, "SLE_pDC_Cleaned_For_pySCENIC.h5ad"))
    print("SUCCESS: Preprocessing Complete.")
