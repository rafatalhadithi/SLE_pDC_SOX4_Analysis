import h5py
import pandas as pd
import scanpy as sc
from arboreto.algo import grnboost2
from arboreto.utils import load_tf_names
from dask.distributed import Client, LocalCluster # <-- NEW: Memory management

if __name__ == '__main__':
    
    print("0. Surgically removing corrupted metadata...")
    with h5py.File("SLE_pDC_Cleaned_For_pySCENIC.h5ad", "r+") as f:
        if 'uns' in f:
            del f['uns']

    print("1. Loading 5,013 pDCs into memory...")
    adata = sc.read_h5ad("SLE_pDC_Cleaned_For_pySCENIC.h5ad")
    ex_matrix = adata.to_df()
    
    print("2. Loading Human Transcription Factor list...")
    tf_names = load_tf_names("hs_hgnc_tfs.txt")
    
    intersecting_tfs = [tf for tf in tf_names if tf in ex_matrix.columns]
    print(f"--> SUCCESS: Found {len(intersecting_tfs)} matching Transcription Factors!")
    
    if len(intersecting_tfs) == 0:
        print("ERROR: 0 matching TFs.")
        exit()
        
    print("3. Spinning up controlled Dask cluster (Max 6 workers, 24GB RAM limit)...")
    # We strictly limit Dask so it cannot crash your 32GB PC
    cluster = LocalCluster(n_workers=6, threads_per_worker=1, memory_limit='4GB')
    custom_client = Client(cluster)
    
    print("4. Starting GRNBoost2 Machine Learning...")
    # We pass the custom_client to arboreto so it obeys our memory rules
    adjacencies = grnboost2(expression_data=ex_matrix, 
                            tf_names=tf_names, 
                            client_or_address=custom_client)
    
    print("5. Calculation complete! Saving results...")
    adjacencies.to_csv("1_adjacencies_output.csv", index=False)
    
    # Safely shut down the memory cluster
    custom_client.close()
    cluster.close()
    print("SUCCESS: Phase 1 is finished.")