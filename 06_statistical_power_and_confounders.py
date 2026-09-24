#!/usr/bin/env python3
import os, scanpy as sc, pandas as pd, numpy as np, scipy.stats as stats, pingouin as pg, matplotlib.pyplot as plt, seaborn as sns

BASE_DIR = os.getcwd()
DATA_DIR, FIGURES_DIR = os.path.join(BASE_DIR, "data"), os.path.join(BASE_DIR, "figures")

if __name__ == "__main__":
    adata = sc.read_h5ad(os.path.join(DATA_DIR, "SLE_pDC_pySCENIC_Complete.h5ad"))
    ifn_genes = [g for g in ['ISG15', 'IFI44', 'IFI44L', 'IFIT1', 'RSAD2', 'IFIT3', 'OAS1'] if g in adata.var_names]
    sc.tl.score_genes(adata, gene_list=ifn_genes, score_name='IFN_Score')

    def get_expr(g): val = adata[:, g].X; return val.toarray().flatten() if not isinstance(val, np.ndarray) else val.flatten()

    df_cells = pd.DataFrame({'Patient': adata.obs['donor_id'], 'SOX4': get_expr('SOX4'), 'CXCR4': get_expr('CXCR4'), 'IFN_Score': adata.obs['IFN_Score']})
    
    pseudobulk_df = df_cells.groupby('Patient').mean().reset_index()
    partial_corr = pg.partial_corr(data=pseudobulk_df, x='SOX4', y='CXCR4', covar='IFN_Score', method='spearman')
    print("PSEUDOBULK PARTIAL CORRELATION:\n", partial_corr)

    np.random.seed(42)
    sox4_full, cxcr4_full = df_cells['SOX4'].values, df_cells['CXCR4'].values
    r_values = []
    
    for _ in range(1000):
        idx = np.random.choice(len(sox4_full), size=127, replace=False)
        corr, pval = stats.spearmanr(sox4_full[idx], cxcr4_full[idx])
        r_values.append(corr)

    plt.figure(figsize=(8, 5))
    sns.histplot(r_values, bins=30, kde=True, color='purple', edgecolor='black')
    plt.axvline(x=0.28, color='red', linestyle='--', linewidth=2, label='Original Effect Size (R=0.28)')
    plt.title('Simulated SOX4-CXCR4 Correlations in Blood (n=127 cells)', fontsize=14, fontweight='bold')
    plt.xlabel('Spearman R Value', fontsize=12); plt.ylabel('Frequency', fontsize=12); plt.legend()
    plt.savefig(os.path.join(FIGURES_DIR, "Figure_5_Bootstrapping_Power_Analysis.pdf"), format='pdf', bbox_inches='tight')
