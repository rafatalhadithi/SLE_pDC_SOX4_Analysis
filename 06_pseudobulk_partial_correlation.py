import pandas as pd
import pingouin as pg
import scanpy as sc
import numpy as np

print("1. Loading SLE pDC dataset...")
adata = sc.read_h5ad("SLE_pDC_pySCENIC_Complete.h5ad")

print("2. Calculating Type I IFN Signature Score on the fly...")
# Canonical gold-standard Type I Interferon signature genes for SLE
ifn_genes = ['ISG15', 'IFI44', 'IFI44L', 'IFIT1', 'RSAD2', 'IFIT3', 'OAS1']
# Ensure genes are in the dataset to avoid errors
ifn_genes_present = [gene for gene in ifn_genes if gene in adata.var_names]
sc.tl.score_genes(adata, gene_list=ifn_genes_present, score_name='IFN_Score')

print("3. Extracting single-cell data...")
# Function to safely extract gene expression whether sparse or dense
def get_expression(adata, gene):
    val = adata[:, gene].X
    return val.toarray().flatten() if not isinstance(val, np.ndarray) else val.flatten()

# Build the dataframe using the correct 'donor_id' column
df = pd.DataFrame({
    'Patient': adata.obs['donor_id'],
    'SOX4': get_expression(adata, 'SOX4'),
    'CXCR4': get_expression(adata, 'CXCR4'),
    'IFN_Score': adata.obs['IFN_Score']
})

print("4. Aggregating into Donor-Level Pseudobulk (Mean expression per patient)...")
pseudobulk_df = df.groupby('Patient').mean().reset_index()
print(f"Aggregated {len(df)} single cells into {len(pseudobulk_df)} donor profiles.")

print("5. Calculating rigorous Pseudobulk Partial Correlation (controlling for IFN storm)...")
partial_corr = pg.partial_corr(data=pseudobulk_df, x='SOX4', y='CXCR4', covar='IFN_Score', method='spearman')

print("\n=======================================================")
print("   PSEUDOBULK PARTIAL CORRELATION RESULTS")
print("=======================================================")
print(partial_corr)
