import scanpy as sc
import liana as li
import pandas as pd

print("Reading metadata to identify cells...")
adata = sc.read_h5ad("SLE_pDC_Cleaned_For_pySCENIC.h5ad")

# Ensure cell types are correctly formatted as strings
adata.obs['cell_type'] = adata.obs['cell_type'].astype(str)

print("Running LIANA Cell-Cell Communication analysis...")
# Run LIANA rank aggregate using consensus resources
li.method.sc.liana_pipe(adata, groupby='cell_type', resource_name='consensus')

# Extract results
results = adata.uns['liana_res']

print("\n==================================================================")
print("  DISCOVERY: CELLS SENDING SIGNALS TO pDC CXCR4 RECEPTORS IN SLE")
print("==================================================================")

# Filter for interactions where pDC is the receiver and CXCR4 is the receptor
pdc_cxcr4_interactions = results[
    (results['target'] == 'pDC') & 
    (results['receptor_complex'] == 'CXCR4') & 
    (results['aggregate_rank'] <= 0.05) # Significant interactions only
]

if pdc_cxcr4_interactions.empty:
    print("No significant CXCR4 interactions found for pDCs.")
else:
    print(pdc_cxcr4_interactions[['source', 'ligand_complex', 'aggregate_rank']])