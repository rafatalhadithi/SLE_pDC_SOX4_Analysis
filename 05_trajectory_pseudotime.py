#!/usr/bin/env python3
import os, scanpy as sc, numpy as np, pandas as pd, matplotlib.pyplot as plt, seaborn as sns
from scipy.stats import spearmanr

BASE_DIR = os.getcwd()
DATA_DIR, FIGURES_DIR = os.path.join(BASE_DIR, "data"), os.path.join(BASE_DIR, "figures")

if __name__ == "__main__":
    adata = sc.read_h5ad(os.path.join(DATA_DIR, "SLE_pDC_pySCENIC_Complete.h5ad"))

    sox4_expr = adata[:, 'SOX4'].X.toarray().flatten() if hasattr(adata.X, 'toarray') else adata[:, 'SOX4'].X.flatten()
    cxcr4_expr = adata[:, 'CXCR4'].X.toarray().flatten() if hasattr(adata.X, 'toarray') else adata[:, 'CXCR4'].X.flatten()
    
    adata.uns['iroot'] = np.argmin(sox4_expr + cxcr4_expr)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.diffmap(adata)
    sc.tl.dpt(adata)

    df = pd.DataFrame({'Pseudotime': adata.obs['dpt_pseudotime'], 'SOX4': sox4_expr, 'CXCR4': cxcr4_expr}).sort_values(by='Pseudotime')
    window = int(len(df) * 0.05) 
    df['SOX4_smoothed'] = df['SOX4'].rolling(window=window, min_periods=1).mean()
    df['CXCR4_smoothed'] = df['CXCR4'].rolling(window=window, min_periods=1).mean()

    plt.figure(figsize=(10, 6))
    sns.lineplot(x='Pseudotime', y='SOX4_smoothed', data=df, label='SOX4 Expression', color='#1f77b4', linewidth=2.5)
    sns.lineplot(x='Pseudotime', y='CXCR4_smoothed', data=df, label='CXCR4 Expression', color='#d62728', linewidth=2.5)
    plt.title('Dynamic Activation of SOX4 and CXCR4 along SLE pDC Trajectory', fontsize=14, fontweight='bold')
    plt.xlabel('Diffusion Pseudotime (Resting → Migratory State)', fontsize=12)
    plt.ylabel('Smoothed Gene Expression', fontsize=12)
    plt.fill_between(df['Pseudotime'], df['SOX4_smoothed'], alpha=0.1, color='#1f77b4')
    plt.fill_between(df['Pseudotime'], df['CXCR4_smoothed'], alpha=0.1, color='#d62728')
    
    plt.savefig(os.path.join(FIGURES_DIR, "Figure_4_Pseudotime_Dynamics.pdf"), format='pdf', bbox_inches='tight')
