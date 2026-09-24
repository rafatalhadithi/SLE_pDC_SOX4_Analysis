#!/usr/bin/env python3
import os, sys, gc, scanpy as sc, liana as li, pandas as pd

BASE_DIR = os.getcwd()
DATA_DIR, OUTPUT_DIR = os.path.join(BASE_DIR, "data"), os.path.join(BASE_DIR, "outputs")

if __name__ == "__main__":
    adata_full = sc.read_h5ad(os.path.join(DATA_DIR, "4118e166-34f5-4c1f-9eed-c64b90a3dace.h5ad"), backed='r')
    
    sle_obs = adata_full.obs[adata_full.obs['disease'] == 'systemic lupus erythematosus']
    pdc_obs = sle_obs[sle_obs['cell_type'] == 'plasmacytoid dendritic cell']
    other_obs = sle_obs[sle_obs['cell_type'] != 'plasmacytoid dendritic cell']
    
    other_sampled = other_obs.sample(n=min(15000, len(other_obs)), random_state=42)
    keep_barcodes = pdc_obs.index.tolist() + other_sampled.index.tolist()
    
    adata_sle = adata_full[keep_barcodes].to_memory()
    del adata_full; gc.collect()

    if 'feature_name' in adata_sle.var.columns: adata_sle.var_names = adata_sle.var['feature_name'].values
    adata_sle.var_names_make_unique()
    if adata_sle.raw is not None: del adata_sle.raw

    li.mt.rank_aggregate(adata_sle, groupby='cell_type', use_raw=False)
    liana_res = adata_sle.uns['liana_res']
    liana_res.to_csv(os.path.join(OUTPUT_DIR, "LIANA_Full_Results.csv"), index=False)
    
    rec_col = 'receptor_complex' if 'receptor_complex' in liana_res.columns else 'receptor'
    lig_col = 'ligand_complex' if 'ligand_complex' in liana_res.columns else 'ligand'
    
    pdc_cxcr4_interactions = liana_res[(liana_res['target'] == 'plasmacytoid dendritic cell') & (liana_res[rec_col].astype(str).str.contains('CXCR4'))]

    if not pdc_cxcr4_interactions.empty:
        top_interactions = pdc_cxcr4_interactions.sort_values(by='magnitude_rank')
        top_interactions.to_csv(os.path.join(OUTPUT_DIR, "pDC_CXCR4_Interactions.csv"), index=False)
    print("SUCCESS: LIANA analysis complete.")
