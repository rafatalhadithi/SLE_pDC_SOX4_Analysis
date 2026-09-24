#!/usr/bin/env python3
import os, glob, pandas as pd, scanpy as sc, scipy.stats as stats, anndata as ad, warnings
warnings.filterwarnings("ignore")

BASE_DIR = os.getcwd()
DATA_DIR = os.path.join(BASE_DIR, "data")

def validate_pediatric_blood():
    print("--- 1. PEDIATRIC SLE BLOOD VALIDATION (GSE135779) ---")
    data_dir = os.path.join(DATA_DIR, "GSE135779_RAW_Data")
    genes_file = os.path.join(DATA_DIR, "GSE135779_genes.tsv.gz")
    
    if not os.path.exists(data_dir) or not os.path.exists(genes_file):
        print("Pediatric data not found. Reviewers with raw matrices can execute this block.")
        return

    genes_df = pd.read_csv(genes_file, sep='\t', header=None)
    gene_symbols = pd.Index(genes_df.iloc[:, 1].values if genes_df.shape[1] > 1 else genes_df.iloc[:, 0].values).astype(str)
    
    adatas = []
    for mtx_path in glob.glob(os.path.join(data_dir, "*_matrix.mtx.gz")):
        adata = sc.read_mtx(mtx_path).T
        adata.obs_names = pd.read_csv(mtx_path.replace("_matrix.mtx.gz", "_barcodes.tsv.gz"), sep='\t', header=None)[0].values
        adata.var_names = gene_symbols
        adata.var_names_make_unique()
        adatas.append(adata)
        
    adata_pediatric = ad.concat(adatas, join="inner")
    sc.pp.normalize_total(adata_pediatric, target_sum=1e4)
    sc.pp.log1p(adata_pediatric)
    
    if 'CLEC4C' in adata_pediatric.var_names and 'LILRA4' in adata_pediatric.var_names:
        c4 = adata_pediatric[:, 'CLEC4C'].X.toarray().flatten() if hasattr(adata_pediatric.X, 'toarray') else adata_pediatric[:, 'CLEC4C'].X.flatten()
        l4 = adata_pediatric[:, 'LILRA4'].X.toarray().flatten() if hasattr(adata_pediatric.X, 'toarray') else adata_pediatric[:, 'LILRA4'].X.flatten()
        
        adata_pdc = adata_pediatric[(c4 > 0) | (l4 > 0)].copy()
        s4 = adata_pdc[:, 'SOX4'].X.toarray().flatten() if hasattr(adata_pdc.X, 'toarray') else adata_pdc[:, 'SOX4'].X.flatten()
        cx4 = adata_pdc[:, 'CXCR4'].X.toarray().flatten() if hasattr(adata_pdc.X, 'toarray') else adata_pdc[:, 'CXCR4'].X.flatten()
        corr, pval = stats.spearmanr(s4, cx4)
        print(f"Isolated {adata_pdc.n_obs} pDCs | R: {corr:.4f}, p: {pval:.2e}")

def validate_adult_kidney():
    print("--- 2. LUPUS NEPHRITIS KIDNEY VALIDATION (GSE303481) ---")
    data_dir = os.path.join(DATA_DIR, "GSE303481_RAW_Data")
    h5_files = glob.glob(os.path.join(data_dir, "*_P*_filtered_feature_bc_matrix.h5"))
    
    if not h5_files:
        print("Kidney matrices not found in data/ folder. Skipping execution.")
        return

    adatas = [sc.read_10x_h5(f) for f in h5_files]
    for a in adatas: a.var_names_make_unique()
    adata_kidney = ad.concat(adatas, join='outer')
    adata_kidney.obs_names_make_unique()
    
    sc.pp.normalize_total(adata_kidney, target_sum=1e4)
    sc.pp.log1p(adata_kidney)
    
    if 'CLEC4C' in adata_kidney.var_names and 'LILRA4' in adata_kidney.var_names:
        c4 = adata_kidney[:, 'CLEC4C'].X.toarray().flatten() if hasattr(adata_kidney.X, 'toarray') else adata_kidney[:, 'CLEC4C'].X.flatten()
        l4 = adata_kidney[:, 'LILRA4'].X.toarray().flatten() if hasattr(adata_kidney.X, 'toarray') else adata_kidney[:, 'LILRA4'].X.flatten()
        
        adata_pdc = adata_kidney[(c4 > 0) | (l4 > 0)].copy()
        s4 = adata_pdc[:, 'SOX4'].X.toarray().flatten() if hasattr(adata_pdc.X, 'toarray') else adata_pdc[:, 'SOX4'].X.flatten()
        cx4 = adata_pdc[:, 'CXCR4'].X.toarray().flatten() if hasattr(adata_pdc.X, 'toarray') else adata_pdc[:, 'CXCR4'].X.flatten()
        corr, pval = stats.spearmanr(s4, cx4)
        print(f"Isolated {adata_pdc.n_obs} pDCs | R: {corr:.4f}, p: {pval:.2e}")

if __name__ == "__main__":
    validate_pediatric_blood()
    validate_adult_kidney()
