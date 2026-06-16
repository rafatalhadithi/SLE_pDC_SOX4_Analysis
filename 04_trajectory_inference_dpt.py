import scanpy as sc
import numpy as np
from scipy.stats import spearmanr

print("1. Loading SLE pDC dataset...")
adata = sc.read_h5ad("SLE_pDC_pySCENIC_Complete.h5ad")

print("2. Calculating robust Naïve/Resting pDC Signature...")
# Canonical resting/steady-state pDC markers
resting_genes = ['TCF4', 'LILRA4', 'CLEC4C', 'IL3RA']
resting_present = [g for g in resting_genes if g in adata.var_names]

# Score cells for this resting signature to smooth out single-gene outliers
sc.tl.score_genes(adata, gene_list=resting_present, score_name='Resting_Score')

# Find the mathematically stable root cell
root_idx = np.argmax(adata.obs['Resting_Score'].values)
print(f"Robust unbiased root cell selected! Index: {root_idx}")

print("3. Running Diffusion Pseudotime (DPT)...")
if 'X_pca' not in adata.obsm:
    sc.tl.pca(adata)
if 'distances' not in adata.obsp:
    sc.pp.neighbors(adata)
if 'X_diffmap' not in adata.obsm:
    sc.tl.diffmap(adata)

adata.uns['iroot'] = root_idx
sc.tl.dpt(adata)

print("4. Evaluating SOX4 and CXCR4 dynamics along the stable trajectory...")
def get_expr(gene):
    expr = adata[:, gene].X
    return expr.toarray().flatten() if not isinstance(expr, np.ndarray) else expr.flatten()

sox4_expr = get_expr('SOX4')
cxcr4_expr = get_expr('CXCR4')
dpt_values = adata.obs['dpt_pseudotime'].values

sox4_r, sox4_p = spearmanr(dpt_values, sox4_expr)
cxcr4_r, cxcr4_p = spearmanr(dpt_values, cxcr4_expr)

print("\n=======================================================")
print("   ROBUST UNBIASED TRAJECTORY RESULTS")
print("=======================================================")
print(f"SOX4 correlation with pseudotime:  R = {sox4_r:.3f}, p = {sox4_p:.2e}")
print(f"CXCR4 correlation with pseudotime: R = {cxcr4_r:.3f}, p = {cxcr4_p:.2e}")
print("=======================================================")
