#!/usr/bin/env python3
"""
Script: 02_pyscenic_grn_inference.py
Purpose: Executes GRNBoost2, RcisTarget, and AUCell to infer regulatory networks, 
         then merges the activity scores into the AnnData object.
"""
import os
import sys
import h5py
import pandas as pd
import scanpy as sc
import ast
from arboreto.algo import grnboost2
from arboreto.utils import load_tf_names
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
from ctxcore.rnkdb import FeatherRankingDatabase
from dask.distributed import Client, LocalCluster

# --- Dynamic Paths ---
BASE_DIR = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
DATA_DIR = os.path.join(BASE_DIR, "data")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")
os.makedirs(OUTPUT_DIR, exist_ok=True)

INPUT_H5AD = os.path.join(DATA_DIR, "SLE_pDC_Cleaned_For_pySCENIC.h5ad")
TF_LIST = os.path.join(DATA_DIR, "hs_hgnc_tfs.txt")
FEATHER_DB = os.path.join(DATA_DIR, "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather")
MOTIF_TBL = os.path.join(DATA_DIR, "motifs-v9-nr.hgnc-m0.001-o0.0.tbl")
FINAL_H5AD = os.path.join(DATA_DIR, "SLE_pDC_pySCENIC_Complete.h5ad")

if __name__ == '__main__':
    print("0. Surgically removing corrupted metadata...")
    with h5py.File(INPUT_H5AD, "r+") as f:
        if 'uns' in f:
            del f['uns']

    print("1. Loading pDCs into memory...")
    adata = sc.read_h5ad(INPUT_H5AD)
    ex_matrix = adata.to_df()
    
    tf_names = load_tf_names(TF_LIST)
    intersecting_tfs = [tf for tf in tf_names if tf in ex_matrix.columns]
    
    print("2. Starting GRNBoost2 Machine Learning...")
    cluster = LocalCluster(n_workers=6, threads_per_worker=1, memory_limit='4GB')
    custom_client = Client(cluster)
    adjacencies = grnboost2(expression_data=ex_matrix, tf_names=intersecting_tfs, client_or_address=custom_client)
    
    print("3. Starting RcisTarget Motif Enrichment...")
    db = FeatherRankingDatabase(fname=FEATHER_DB, name="hg38")
    motifs = prune2df(rnkdbs=[db], modules=list(adjacencies), motif_annotations_fname=MOTIF_TBL, client_or_address=custom_client)
    motifs.to_csv(os.path.join(OUTPUT_DIR, "Table_S1_Regulons_Raw.csv"))
    
    print("4. Formatting Regulons & Running AUCell...")
    # Fix pandas multi-index quirk
    motifs[('Enrichment', 'TargetGenes')] = motifs[('Enrichment', 'TargetGenes')].apply(lambda x: eval(x) if isinstance(x, str) else x)
    regulons = df2regulons(motifs)
    auc_mtx = aucell(ex_matrix, regulons, num_workers=4)
    
    print("5. Merging & Saving Final Dataset...")
    adata.obs = adata.obs.join(auc_mtx)
    adata.write_h5ad(FINAL_H5AD)
    
    print("6. Extracting IRF7 Target Genes for GO Enrichment...")
    for reg in regulons:
        if reg.name == 'IRF7(+)':
            with open(os.path.join(DATA_DIR, "IRF7_Target_Genes.txt"), "w") as f:
                f.write("\n".join(list(reg.genes)))
                
    custom_client.close()
    cluster.close()
    print("SUCCESS: pySCENIC Pipeline Complete!")
