#!/usr/bin/env python3
import os, scanpy as sc, harmonypy as hm

BASE_DIR = os.getcwd()
DATA_DIR = os.path.join(BASE_DIR, "data")

if __name__ == "__main__":
    adata = sc.read_h5ad(os.path.join(DATA_DIR, "SLE_pDC_pySCENIC_Complete.h5ad"))
    sc.tl.pca(adata, svd_solver='arpack')

    ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, ['donor_id'])
    adata.obsm['X_pca_harmony'] = ho.Z_corr if ho.Z_corr.shape[0] == adata.n_obs else ho.Z_corr.T

    sc.pp.neighbors(adata, use_rep='X_pca_harmony')
    sc.tl.umap(adata)
    adata.write_h5ad(os.path.join(DATA_DIR, "SLE_pDC_Integrated.h5ad"))
    print("SUCCESS: Harmony integration complete.")
