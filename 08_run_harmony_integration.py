import scanpy as sc
import matplotlib.pyplot as plt
import harmonypy as hm

print("1. Loading Adult pDC dataset...")
adata = sc.read_h5ad("SLE_pDC_pySCENIC_Complete.h5ad")

print(f"Total cells loaded: {adata.n_obs}")
print(f"Total distinct donors to integrate: {adata.obs['donor_id'].nunique()}")

print("2. Computing standard PCA...")
sc.tl.pca(adata, svd_solver='arpack')

print("3. Running native Harmony integration...")
# Call harmonypy directly to bypass Scanpy wrapper bugs
ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, ['donor_id'])

# BULLETPROOF MATRICE ASSIGNMENT:
# Automatically detect matrix orientation and assign safely
harmony_matrix = ho.Z_corr
if harmony_matrix.shape[0] == adata.n_obs:
    adata.obsm['X_pca_harmony'] = harmony_matrix
else:
    adata.obsm['X_pca_harmony'] = harmony_matrix.T

print("4. Re-calculating nearest neighbors and UMAP on the newly integrated data...")
sc.pp.neighbors(adata, use_rep='X_pca_harmony')
sc.tl.umap(adata)

print("5. Plotting the Integrated UMAP...")
sc.pl.umap(adata, color=['donor_id', 'disease', 'SOX4'], ncols=3, show=False)
plt.savefig("Harmony_Integrated_UMAP.png", bbox_inches='tight', dpi=300)

print("6. Saving Harmony-integrated dataset...")
adata.write_h5ad("SLE_pDC_Integrated.h5ad")

print("\n>>> MASSIVE SUCCESS! Data is integrated. <<<")
print("Saved plot as: Harmony_Integrated_UMAP.png")
print("Saved data as: SLE_pDC_Integrated.h5ad")
