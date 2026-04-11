import scanpy as sc
import pandas as pd
import scipy.stats as stats
import pingouin as pg

print("Loading annotated SLE pDC dataset...")
adata = sc.read_h5ad("SLE_pDC_Cleaned_For_pySCENIC.h5ad")

# --- 1. DONOR-LEVEL PSEUDOBULK CORRELATION ---
print("Aggregating cells to donor-level pseudobulks...")
# Group by patient donor ID and calculate mean expression
pseudobulk = adata.to_df().groupby(adata.obs['donor_id']).mean()

sox4_bulk = pseudobulk['SOX4']
cxcr4_bulk = pseudobulk['CXCR4']

corr_pb, pval_pb = stats.spearmanr(sox4_bulk, cxcr4_bulk)
print(f"Pseudobulk Spearman R: {corr_pb:.2f}, p-value: {pval_pb:.2e}")

# --- 2. PARTIAL CORRELATION (REGRESSING OUT IFN SIGNATURE) ---
print("Calculating partial correlation controlling for IFN score...")
# Assuming an 'IFN_Score' was previously calculated and stored in obs
df_stats = pd.DataFrame({
    'SOX4': adata[:, 'SOX4'].X.toarray().flatten(),
    'CXCR4': adata[:, 'CXCR4'].X.toarray().flatten(),
    'IFN_Score': adata.obs['IFN_Score'].values
})

# Use pingouin to run partial spearman
partial_corr = pg.partial_corr(data=df_stats, x='SOX4', y='CXCR4', covar='IFN_Score', method='spearman')
print(partial_corr)