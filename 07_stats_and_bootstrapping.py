#!/usr/bin/env python3
"""
Script: 07_statistical_power_and_confounders.py
Purpose: Evaluates pseudobulk partial correlations controlling for IFN signatures 
         and models Type II error rates via bootstrapping simulation.
"""
import os
import sys
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.stats as stats
import pingouin as pg
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings("ignore")

BASE_DIR = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
DATA_DIR = os.path.join(BASE_DIR, "data")
FIGURES_DIR = os.path.join(BASE_DIR, "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)

INPUT_H5AD = os.path.join(DATA_DIR, "SLE_pDC_pySCENIC_Complete.h5ad")

if __name__ == "__main__":
    if not os.path.exists(INPUT_H5AD):
        print(f"❌ ERROR: {INPUT_H5AD} not found.")
        sys.exit(1)

    print("1. Loading annotated SLE pDC dataset...")
    adata = sc.read_h5ad(INPUT_H5AD)

    print("2. Calculating Type I IFN Signature Score...")
    ifn_genes = ['ISG15', 'IFI44', 'IFI44L', 'IFIT1', 'RSAD2', 'IFIT3', 'OAS1']
    ifn_genes_present = [g for g in ifn_genes if g in adata.var_names]
    sc.tl.score_genes(adata, gene_list=ifn_genes_present, score_name='IFN_Score')

    def get_expression(adata, gene):
        val = adata[:, gene].X
        return val.toarray().flatten() if not isinstance(val, np.ndarray) else val.flatten()

    df_cells = pd.DataFrame({
        'Patient': adata.obs['donor_id'],
        'SOX4': get_expression(adata, 'SOX4'),
        'CXCR4': get_expression(adata, 'CXCR4'),
        'IFN_Score': adata.obs['IFN_Score']
    })

    print("3. Aggregating into Donor-Level Pseudobulk (n=242)...")
    pseudobulk_df = df_cells.groupby('Patient').mean().reset_index()

    print("4. Calculating Pseudobulk Partial Correlation (controlling for IFN storm)...")
    partial_corr = pg.partial_corr(data=pseudobulk_df, x='SOX4', y='CXCR4', covar='IFN_Score', method='spearman')
    print("\n--- PSEUDOBULK PARTIAL CORRELATION RESULTS ---")
    print(partial_corr.to_string())

    print("\n5. Running Bootstrapping Power Simulation (1,000 iterations at n=127)...")
    np.random.seed(42)
    sox4_full = df_cells['SOX4'].values
    cxcr4_full = df_cells['CXCR4'].values
    
    significant_count = 0
    r_values = []
    
    for _ in range(1000):
        random_indices = np.random.choice(len(sox4_full), size=127, replace=False)
        corr, pval = stats.spearmanr(sox4_full[random_indices], cxcr4_full[random_indices])
        r_values.append(corr)
        if pval < 0.05 and corr > 0:
            significant_count += 1

    power_percentage = (significant_count / 1000) * 100
    print(f"Statistical Power to detect the network: {power_percentage:.1f}%")

    plt.figure(figsize=(8, 5))
    plt.rcParams.update({'font.sans-serif': 'Arial'})
    sns.histplot(r_values, bins=30, kde=True, color='purple', edgecolor='black')
    plt.axvline(x=0.28, color='red', linestyle='--', linewidth=2, label='Original Blood Effect Size (R=0.28)')
    plt.title('Simulated SOX4-CXCR4 Correlations in Blood (n=127 cells)', fontsize=14, fontweight='bold')
    plt.xlabel('Spearman R Value', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.legend()
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    plt.savefig(os.path.join(FIGURES_DIR, "Figure_5_Bootstrapping_Power_Analysis.pdf"), format='pdf', bbox_inches='tight')
    print("SUCCESS: Simulation complete and plotted.")
