#!/usr/bin/env python3
"""
Script: 06_liana_cell_communication.py
Purpose: Uses the LIANA framework to map ligand-receptor interactions across the PBMC cohort, 
         specifically screening for peripheral CXCR4 targeting.
"""
import os
import sys
import gc
import scanpy as sc
import liana as li
import pandas as pd

BASE_DIR = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
DATA_DIR = os.path.join(BASE_DIR, "data")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")
os.makedirs(OUTPUT_DIR, exist_ok=True)

GLOBAL_H5AD = os.path.join(DATA_DIR, "4118e166-34f5-4c1f-9eed-c64b90a3dace.h5ad")

if __name__ == "__main__":
    if not os.path.exists(GLOBAL_H5AD):
        print(f"❌ ERROR: Global atlas not found at {GLOBAL_H5AD}")
        sys.exit(1)

    print("1. Reading metadata to isolate SLE cells...")
    adata_full = sc.read_h5ad(GLOBAL_H5AD, backed='r')
    obs_df = adata_full.obs
    
    sle_obs = obs_df[obs_df['disease'] == 'systemic lupus erythematosus']
    pdc_obs = sle_obs[sle_obs['cell_type'] == 'plasmacytoid dendritic cell']
    other_obs = sle_obs[sle_obs['cell_type'] != 'plasmacytoid dendritic cell']
    
    print("2. Subsampling peripheral immune cells to conserve memory...")
    n_sample = min(15000, len(other_obs))
    other_sampled = other_obs.sample(n=n_sample, random_state=42)
    keep_barcodes = pdc_obs.index.tolist() + other_sampled.index.tolist()
    
    adata_sle = adata_full[keep_barcodes].to_memory()
    del adata_full
    gc.collect()

    print("3. Converting Ensembl IDs to Gene Symbols...")
    if 'feature_name' in adata_sle.var.columns:
        adata_sle.var_names = adata_sle.var['feature_name'].values
    adata_sle.var_names_make_unique()
    
    if adata_sle.raw is not None:
        del adata_sle.raw

    print("4. Running LIANA Cell-Cell Communication analysis...")
    li.mt.rank_aggregate(adata_sle, groupby='cell_type', use_raw=False)
    
    liana_res = adata_sle.uns['liana_res']
    liana_res.to_csv(os.path.join(OUTPUT_DIR, "LIANA_Full_Results.csv"), index=False)
    
    rec_col = 'receptor_complex' if 'receptor_complex' in liana_res.columns else 'receptor'
    lig_col = 'ligand_complex' if 'ligand_complex' in liana_res.columns else 'ligand'
    
    pdc_cxcr4_interactions = liana_res[
        (liana_res['target'] == 'plasmacytoid dendritic cell') & 
        (liana_res[rec_col].astype(str).str.contains('CXCR4'))
    ]

    print("\n--- DISCOVERY: CELLS SENDING SIGNALS TO pDC CXCR4 RECEPTORS IN SLE ---")
    if pdc_cxcr4_interactions.empty:
        print("Result: No significant CXCR4 interactions found for pDCs in systemic circulation.")
    else:
        top_interactions = pdc_cxcr4_interactions.sort_values(by='magnitude_rank')
        print(top_interactions[['source', lig_col, 'target', rec_col, 'magnitude_rank']].head(10).to_string(index=False))
        top_interactions.to_csv(os.path.join(OUTPUT_DIR, "pDC_CXCR4_Interactions.csv"), index=False)
    print("SUCCESS: LIANA analysis complete.")
