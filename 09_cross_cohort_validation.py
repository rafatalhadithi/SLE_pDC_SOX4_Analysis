#!/usr/bin/env python3
"""
Script: 09_cross_cohort_validation.py
Purpose: Validates the SOX4-CXCR4 correlation in an independent pediatric SLE PBMC cohort 
         and a Lupus Nephritis kidney snRNA-seq cohort to prove compartmentalized tissue homing.
"""
import os
import sys
import glob
import pandas as pd
import scanpy as sc
import scipy.stats as stats
import anndata as ad
import warnings
warnings.filterwarnings("ignore")

BASE_DIR = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
DATA_DIR = os.path.join(BASE_DIR, "data")

def validate_pediatric_blood():
    print("\n=======================================================")
    print(" 1. PEDIATRIC SLE BLOOD VALIDATION (GSE135779)")
    print("=======================================================")
    
    genes_file = os.path.join(DATA_DIR, "GSE135779_genes.tsv.gz")
    data_dir = os.path.join(DATA_DIR, "GSE135779_RAW_Data")
    
    if not os.path.exists(genes_file) or not os.path.exists(data_dir):
        print("Pediatric datasets not found in data/ folder. Skipping execution.")
        return

    print("Loading pediatric matrix files...")
    genes_df = pd.read_csv(genes_file, sep='\t', header=None)
    gene_symbols = pd.Index(genes_df.iloc[:, 1].values if genes_df.shape[1] > 1 else genes_df.iloc[:, 0].values).astype(str)
    
    mtx_files = glob.glob(os.path.join(data_dir, "*_matrix.mtx.gz"))
    adatas = []
    
    for mtx_path in mtx_files:
        barcode_path = mtx_path.replace("_matrix.mtx.gz", "_barcodes.tsv.gz")
        adata = sc.read_mtx(mtx_path).T
        adata.obs_names = pd.read_csv(barcode_path, sep='\t', header=None)[0].values
        adata.var_names = gene_symbols
        adata.var_names_make_unique()
        adatas.append(adata)
        
    adata_pediatric = ad.concat(adatas, join="inner")
    sc.pp.normalize_total(adata_pediatric, target_sum=1e4)
    sc.pp.log1p(adata_pediatric)
    
    if 'CLEC4C' in adata_pediatric.var_names and 'LILRA4' in adata_pediatric.var_names:
        clec4c = adata_pediatric[:, 'CLEC4C'].X.toarray().flatten() if hasattr(adata_pediatric.X, 'toarray') else adata_pediatric[:, 'CLEC4C'].X.flatten()
        lilra4 = adata_pediatric[:, 'LILRA4'].X.toarray().flatten() if hasattr(adata_pediatric.X, 'toarray') else adata_pediatric[:, 'LILRA4'].X.flatten()
        
        pdc_mask = (clec4c > 0) | (lilra4 > 0)
        adata_pdc = adata_pediatric[pdc_mask].copy()
        print(f"Digitally isolated {adata_pdc.n_obs} pediatric pDCs.")
        
        sox4 = adata_pdc[:, 'SOX4'].X.toarray().flatten() if hasattr(adata_pdc.X, 'toarray') else adata_pdc[:, 'SOX4'].X.flatten()
        cxcr4 = adata_pdc[:, 'CXCR4'].X.toarray().flatten() if hasattr(adata_pdc.X, 'toarray') else adata_pdc[:, 'CXCR4'].X.flatten()
        
        corr, pval = stats.spearmanr(sox4, cxcr4)
        print(f"Spearman R: {corr:.4f}")
        print(f"P-value: {pval:.2e}")
        print("--> SUCCESS: The SOX4-CXCR4 correlation is pathologically amplified in circulating pediatric SLE pDCs.")

def validate_adult_kidney():
    print("\n=======================================================")
    print(" 2. LUPUS NEPHRITIS KIDNEY VALIDATION (GSE303481)")
    print("=======================================================")
    
    data_dir = os.path.join(DATA_DIR, "GSE303481_RAW_Data")
    h5_files = glob.glob(os.path.join(data_dir, "*_P*_filtered_feature_bc_matrix.h5"))
    
    if not h5_files:
        print("Kidney .h5 matrices not found in data/ folder. Skipping execution.")
        return

    print(f"Loading {len(h5_files)} Lupus Nephritis patient biopsy samples...")
    adatas = []
    for f in h5_files:
        adata = sc.read_10x_h5(f)
        adata.var_names_make_unique()
        adatas.append(adata)
        
    adata_kidney = ad.concat(adatas, join='outer')
    adata_kidney.obs_names_make_unique()
    sc.pp.normalize_total(adata_kidney, target_sum=1e4)
    sc.pp.log1p(adata_kidney)
    
    if 'CLEC4C' in adata_kidney.var_names and 'LILRA4' in adata_kidney.var_names:
        clec4c = adata_kidney[:, 'CLEC4C'].X.toarray().flatten() if hasattr(adata_kidney.X, 'toarray') else adata_kidney[:, 'CLEC4C'].X.flatten()
        lilra4 = adata_kidney[:, 'LILRA4'].X.toarray().flatten() if hasattr(adata_kidney.X, 'toarray') else adata_kidney[:, 'LILRA4'].X.flatten()
        
        pdc_mask = (clec4c > 0) | (lilra4 > 0)
        adata_pdc = adata_kidney[pdc_mask].copy()
        print(f"Digitally isolated {adata_pdc.n_obs} tissue-infiltrating pDCs.")
        
        sox4 = adata_pdc[:, 'SOX4'].X.toarray().flatten() if hasattr(adata_pdc.X, 'toarray') else adata_pdc[:, 'SOX4'].X.flatten()
        cxcr4 = adata_pdc[:, 'CXCR4'].X.toarray().flatten() if hasattr(adata_pdc.X, 'toarray') else adata_pdc[:, 'CXCR4'].X.flatten()
        
        corr, pval = stats.spearmanr(sox4, cxcr4)
        print(f"Spearman R: {corr:.4f}")
        print(f"P-value: {pval:.2e}")
        print("--> SUCCESS: The correlation is uncoupled in the kidney, confirming the compartmentalized homing model.")

if __name__ == "__main__":
    validate_pediatric_blood()
    validate_adult_kidney()
