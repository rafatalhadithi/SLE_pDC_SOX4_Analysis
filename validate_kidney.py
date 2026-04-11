import scanpy as sc
import pandas as pd
import glob
import os
import sys
import scipy.stats as stats
import anndata as ad

print("1. Loading universal gene list...")
genes_file = "GSE135779_genes.tsv.gz"
try:
    # Load the master gene file
    genes_df = pd.read_csv(genes_file, sep='\t', header=None)
    # GEO files usually have Ensembl IDs in col 0 and Symbols in col 1
    gene_symbols = genes_df.iloc[:, 1].values if genes_df.shape[1] > 1 else genes_df.iloc[:, 0].values
    gene_symbols = pd.Index(gene_symbols).astype(str)
except Exception as e:
    print(f"Error loading genes: {e}")
    sys.exit(1)

print("2. Finding kidney matrix files...")
data_dir = "GSE135779_RAW_Data"
mtx_files = glob.glob(os.path.join(data_dir, "*_matrix.mtx.gz"))

if not mtx_files:
    print("No matrix files found!")
    sys.exit(1)

print(f"Found {len(mtx_files)} patient samples. Loading data (this takes ~60 seconds)...")
adatas = []
for mtx_path in mtx_files:
    barcode_path = mtx_path.replace("_matrix.mtx.gz", "_barcodes.tsv.gz")
    
    # Load matrix and transpose (genes to columns, cells to rows)
    adata = sc.read_mtx(mtx_path).T
    
    # Attach the barcodes and genes
    barcodes = pd.read_csv(barcode_path, sep='\t', header=None)[0].values
    adata.obs_names = barcodes
    adata.var_names = gene_symbols
    adata.var_names_make_unique()
    
    adatas.append(adata)

print("3. Merging all Lupus Nephritis kidney samples...")
adata_kidney = ad.concat(adatas, join="inner")
print(f"Total kidney cells loaded: {adata_kidney.n_obs}")

print("4. Normalizing and log-transforming expression...")
sc.pp.normalize_total(adata_kidney, target_sum=1e4)
sc.pp.log1p(adata_kidney)

print("5. Isolating pDCs from the kidney infiltrates...")
# CLEC4C (BDCA-2) and LILRA4 are the gold standard human pDC markers
if 'CLEC4C' in adata_kidney.var_names and 'LILRA4' in adata_kidney.var_names:
    clec4c_expr = adata_kidney[:, 'CLEC4C'].X.toarray().flatten()
    lilra4_expr = adata_kidney[:, 'LILRA4'].X.toarray().flatten()
    
    # Keep cells that express either of the pDC markers
    pdc_mask = (clec4c_expr > 0) | (lilra4_expr > 0)
    adata_pdc = adata_kidney[pdc_mask]
    
    print(f"Found {adata_pdc.n_obs} pDCs that successfully migrated into the kidneys.")
    
    if adata_pdc.n_obs > 10:
        sox4_expr = adata_pdc[:, 'SOX4'].X.toarray().flatten()
        cxcr4_expr = adata_pdc[:, 'CXCR4'].X.toarray().flatten()
        
        corr, pval = stats.spearmanr(sox4_expr, cxcr4_expr)
        
        print("\n==========================================================")
        print("  LUPUS NEPHRITIS KIDNEY VALIDATION (SOX4 vs CXCR4)")
        print("==========================================================")
        print(f"Spearman Correlation (R): {corr:.4f}")
        print(f"P-value:                  {pval:.2e}")
        
        if corr > 0 and pval < 0.05:
            print("\n>>> SUCCESS! The SOX4-CXCR4 network is active inside the target organ! <<<")
    else:
        print("Not enough pDCs found to perform correlation.")
else:
    print("Could not find pDC marker genes in this dataset.")