#!/usr/bin/env python3
"""
Script: 08_harmony_batch_integration.py
Purpose: Executes native harmonypy batch-effect correction across heterogeneous donor profiles.
"""
import os
import sys
import scanpy as sc
import harmonypy as hm

BASE_DIR = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
DATA_DIR = os.path.join(BASE_DIR, "data")
INPUT_H5AD = os.path.join(DATA_DIR, "SLE_pDC_pySCENIC_Complete.h5ad")
OUTPUT_H5AD = os.path.join(DATA_DIR, "SLE_pDC_Integrated.h5ad")

if __name__ == "__main__":
    if not os.path.exists(INPUT_H5AD):
        print(f"❌ ERROR: {INPUT_H5AD} not found.")
        sys.exit(1)

    print("1. Loading Adult pDC dataset...")
    adata = sc.read_h5ad(INPUT_H5AD)
    
    print("2. Computing standard PCA...")
    sc.tl.pca(adata, svd_solver='arpack')

    print("3. Running native Harmony integration across donors...")
    ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, ['donor_id'])
    
    harmony_matrix = ho.Z_corr
    if harmony_matrix.shape[0] == adata.n_obs:
        adata.obsm['X_pca_harmony'] = harmony_matrix
    else:
        adata.obsm['X_pca_harmony'] = harmony_matrix.T

    print("4. Re-calculating nearest neighbors and UMAP...")
    sc.pp.neighbors(adata, use_rep='X_pca_harmony')
    sc.tl.umap(adata)

    print("5. Saving Harmony-integrated dataset...")
    adata.write_h5ad(OUTPUT_H5AD)
    print("SUCCESS: Data integrated.")
