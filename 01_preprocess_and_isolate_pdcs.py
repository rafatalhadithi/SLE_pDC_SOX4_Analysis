# 01_preprocess_and_isolate_pdcs.py
# Auto-generated from original preprocessing notebook

import scanpy as sc

# 1. Load the dataset in "backed" mode (using relative path for anonymity and portability)
adata = sc.read_h5ad("data/4118e166-34f5-4c1f-9eed-c64b90a3dace.h5ad", backed='r')

# 2. Filter for Female pDCs and load ONLY them into memory
adata_pdc = adata[(adata.obs['sex'] == 'female') & 
                  (adata.obs['cell_type'] == 'plasmacytoid dendritic cell')].to_memory()

print(f"Successfully loaded {adata_pdc.n_obs} pDCs into memory!")

# 3. Define our XCI-Escape Signature
target_genes = ['XIST', 'KDM6A', 'CXorf38', 'TLR7']

# Verify which genes are present in the matrix
genes_present = [gene for gene in target_genes if gene in adata_pdc.var_names]
print(f"Genes ready for plotting: {genes_present}")

# 4. Generate the Viability Plot (SLE vs Healthy)
if genes_present:
    sc.pl.violin(adata_pdc, keys=genes_present, groupby='disease', rotation=45)
else:
    print("Error: Target genes not found in var_names. We may need to map Ensembl IDs.")

# 1. Print the names of the gene metadata columns
print("--- Gene Metadata Columns ---")
print(adata_pdc.var.columns.tolist())

# 2. Print the first 3 rows to see the exact formatting
print("\n--- First 3 Rows ---")
print(adata_pdc.var.head(3))

import scanpy as sc

# 1. Swap the index from Ensembl IDs to readable Gene Symbols
adata_pdc.var_names = adata_pdc.var['feature_name']

# 2. Make sure all gene names are unique (a strict requirement for Scanpy)
adata_pdc.var_names_make_unique()

# 3. Define our XCI-Escape target genes
target_genes = ['XIST', 'KDM6A', 'CXorf38', 'TLR7']

# 4. Check if they are found now
genes_present = [gene for gene in target_genes if gene in adata_pdc.var_names]
print(f"Genes ready for plotting: {genes_present}")

# 5. Draw the Viability Plot (SLE vs Normal)
if genes_present:
    sc.pl.violin(adata_pdc, keys=genes_present, groupby='disease', rotation=45)
else:
    print("Error: Genes still not found. They may have been filtered out during QC.")

# Force Scanpy to use the main matrix (where we fixed the names) instead of the raw backup
sc.pl.violin(adata_pdc, keys=genes_present, groupby='disease', rotation=45, use_raw=False)

import scanpy as sc

# 1. Exact Ensembl IDs for our human target genes
ensembl_dict = {
    'XIST': 'ENSG00000229807',
    'KDM6A': 'ENSG00000147050',
    'CXorf38': 'ENSG00000183753',
    'TLR7': 'ENSG00000196664'
}

# 2. Find which ones exist in the hidden RAW matrix
raw_keys_to_plot = [eid for name, eid in ensembl_dict.items() if eid in adata_pdc.raw.var_names]

# 3. Print a legend so you know which ID is which on the plot
print("--- PLOT LEGEND ---")
for name, eid in ensembl_dict.items():
    if eid in raw_keys_to_plot:
        print(f"{eid} = {name}")

# 4. Plot directly from the unscaled RAW matrix (use_raw=True)
if raw_keys_to_plot:
    sc.pl.violin(adata_pdc, keys=raw_keys_to_plot, groupby='disease', rotation=45, use_raw=True)
else:
    print("CRITICAL DEAD END: The authors completely scrubbed these genes from the raw data.")

# 1. Normalize and Log-transform the raw data for statistical testing
sc.pp.normalize_total(adata_pdc, target_sum=1e4)
sc.pp.log1p(adata_pdc)

# 2. Run the Differential Gene Expression test (SLE vs Normal)
sc.tl.rank_genes_groups(adata_pdc, groupby='disease', method='wilcoxon', 
                        reference='normal', use_raw=False)

# 3. Extract and print the top 20 significantly upregulated genes in SLE
top_genes = sc.get.rank_genes_groups_df(adata_pdc, group='systemic lupus erythematosus')

# Filter for strong positive fold-change and significance
upregulated_escape_candidates = top_genes[(top_genes['logfoldchanges'] > 0.5) & 
                                          (top_genes['pvals_adj'] < 0.05)]

print("--- TOP NOVEL ESCAPE CANDIDATES IN SLE pDCs ---")
print(upregulated_escape_candidates.head(20)[['names', 'logfoldchanges', 'pvals_adj']])

import scanpy as sc

# 1. Promote the raw, unscaled backup matrix to be our main working object
adata_clean = adata_pdc.raw.to_adata()

# 2. Filter out any "dead" cells with zero total counts to fix the first warning
sc.pp.filter_cells(adata_clean, min_counts=1)

# 3. Safely normalize and log-transform the pristine positive counts
sc.pp.normalize_total(adata_clean, target_sum=1e4)
sc.pp.log1p(adata_clean)

# 4. Run the DEG test (SLE vs Normal)
sc.tl.rank_genes_groups(adata_clean, groupby='disease', method='wilcoxon', reference='normal')

# 5. Extract and print the top significantly upregulated genes in SLE
top_genes = sc.get.rank_genes_groups_df(adata_clean, group='systemic lupus erythematosus')

# Filter for strong positive fold-change and statistical significance
upregulated_escape_candidates = top_genes[(top_genes['logfoldchanges'] > 0.5) & 
                                          (top_genes['pvals_adj'] < 0.05)]

print(f"Total significant upregulated genes found: {len(upregulated_escape_candidates)}")
print("\n--- TOP NOVEL ESCAPE CANDIDATES IN SLE pDCs ---")
print(upregulated_escape_candidates.head(20)[['names', 'logfoldchanges', 'pvals_adj']])

import pandas as pd

# 1. Create a dictionary mapping Ensembl IDs to Gene Symbols from the metadata
gene_map = adata_pdc.var['feature_name'].to_dict()

# 2. Add the Gene Symbols to our results table (using .copy() to prevent pandas warnings)
final_candidates = upregulated_escape_candidates.copy()
final_candidates['Gene_Symbol'] = final_candidates['names'].map(gene_map)

# 3. Reorganize the columns so it's easy to read
final_candidates = final_candidates[['names', 'Gene_Symbol', 'logfoldchanges', 'pvals_adj']]

# 4. Print the translated top 20 genes
print("--- TOP UP-REGULATED GENES IN SLE pDCs (TRANSLATED) ---")
print(final_candidates.head(20))

# Optional: Save the full list of 77 genes to a CSV on your computer for your paper!
final_candidates.to_csv("SLE_pDC_Upregulated_Genes.csv", index=False)

# 1. Build the dictionary from the clean raw object (which still has Ensembl keys)
correct_gene_map = adata_clean.var['feature_name'].to_dict()

# 2. Re-map the symbols using the correct dictionary
final_candidates['Gene_Symbol'] = final_candidates['names'].map(correct_gene_map)

# 3. Print the highly anticipated translated list!
print("--- TOP 20 UP-REGULATED GENES IN SLE pDCs ---")
print(final_candidates.head(20))

# 1. Ensure the var_names are set to the readable Gene Symbols for pySCENIC
adata_clean.var_names = adata_clean.var['feature_name']
adata_clean.var_names_make_unique()

# 2. Save the fully prepped pDC matrix to your hard drive
adata_clean.write_h5ad("SLE_pDC_Cleaned_For_pySCENIC.h5ad")

# 3. Save the full list of 77 upregulated genes for your paper
final_candidates.to_csv("SLE_pDC_Top77_Upregulated.csv", index=False)

print("Export complete! Ready for pySCENIC.")

# 1. Clear the conflicting index name so AnnData allows the save
adata_clean.var.index.name = None

# 2. Save the fully prepped pDC matrix to your hard drive
adata_clean.write_h5ad("SLE_pDC_Cleaned_For_pySCENIC.h5ad")

# 3. Save the full list of 77 upregulated genes for your paper
final_candidates.to_csv("SLE_pDC_Top77_Upregulated.csv", index=False)

print("Export complete! Ready for pySCENIC.")
