# Single-cell Analysis Workflow for SLE pDC Regulon Inference

This repository contains the computational pipeline used to process, integrate, and perform network inference on plasmacytoid dendritic cell (pDC) scRNA-seq data from patients with Systemic Lupus Erythematosus (SLE).

## Overview
This workflow utilizes `Scanpy`, `pySCENIC`, `LIANA`, and `Harmony` to identify the **SOX4-CXCR4** regulatory axis. The analysis tracks the transition of pDCs from a steady state to a hyperactivated, tissue-homing state characterized by infiltration in Lupus Nephritis.

---

## Repository Contents & File I/O Mapping
The scripts are designed to be executed sequentially. Below is the precise input and expected output tracking mapping for each script:

### `01_preprocess_and_isolate_pdcs.py`
Initial quality control, filtering, and isolation of target cell populations.
* **Input File:** `data/4118e166-34f5-4c1f-9eed-c64b90a3dace.h5ad` (Raw large-cohort file)
* **Output Files:** 
  * `SLE_pDC_Cleaned_For_pySCENIC.h5ad` (Isolated high-quality female pDC matrix)
  * `SLE_pDC_Top77_Upregulated.csv` / `SLE_pDC_Upregulated_Genes.csv` (Differentially expressed gene lists)

### `02_run_pyscenic_network.py`
Executes machine learning gene regulatory network inference via GRNBoost2.
* **Input Files:** 
  * `SLE_pDC_Cleaned_For_pySCENIC.h5ad`
  * `hs_hgnc_tfs.txt` (Curated human transcription factor list)
* **Output File:** `1_adjacencies_output.csv` (Inferred regulatory co-expression links)

### `03_go_functional_enrichment.py`
Performs Gene Ontology pathway enrichment testing for the major regulatory modules.
* **Input File:** `IRF7_Target_Genes.txt` (Downstream targets identified inside the network)
* **Output Artifacts:** 
  * `5_IRF7_GO_Enrichment.png` (Pathway bar plot showing significance distribution)
  * Console printout mapping the top 10 significantly enriched biological processes

### `04_trajectory_inference_dpt.py`
Constructs single-cell developmental and activation timelines using Diffusion Pseudotime.
* **Input File:** `SLE_pDC_pySCENIC_Complete.h5ad` (Matrix containing pySCENIC regulon activity metrics)
* **Output Artifacts:** Console metrics displaying Spearman correlation coefficient values (R) and statistical significance (p-values) tracking expression changes along the timeline.

### `05_liana_cell_communication.py`
Maps ligand-receptor interaction networks signaling into the target cell types.
* **Input File:** `SLE_pDC_Cleaned_For_pySCENIC.h5ad`
* **Output Artifacts:** Console printout showing filtered significant upstream source cell types actively signaling to the pDC CXCR4 surface receptor.

### `06_pseudobulk_partial_correlation.py`
Executes donor-level statistical modeling controlling for biological covariates like the systemic interferon storm.
* **Input File:** `SLE_pDC_pySCENIC_Complete.h5ad`
* **Output Artifacts:** A rigorous partial correlation statistics matrix generated at patient-level resolution.

### `07_bootstrapping_power_analysis.py`
Runs a 1,000-iteration down-sampling simulation to evaluate statistical power and control for scRNA-seq technical dropout artifacts.
* **Input File:** `SLE_pDC_Cleaned_For_pySCENIC.h5ad`
* **Output Files:** `Bootstrapping_Power_Analysis.png` (Frequency distribution histogram of correlation boundaries)

### `08_run_harmony_integration.py`
Executes native batch-effect correction across heterogeneous donor background profiles.
* **Input File:** `SLE_pDC_pySCENIC_Complete.h5ad`
* **Output Files:** 
  * `SLE_pDC_Integrated.h5ad` (The standardized, batch-corrected unified multi-donor dataset)
  * `Harmony_Integrated_UMAP.png` (Diagnostic visual control plot showing cross-donor blending)

### `09_plot_publication_umap.py`
Generates high-resolution publication-grade visual figures summarizing global cell attributes.
* **Input File:** `SLE_pDC_Integrated.h5ad`
* **Output File:** `Figure_Harmony_UMAP_Publication.pdf` (High-resolution vector graphic composite showing integration patterns, clinical disease state split, and structural gene expressions side by side)

---

## Data Availability
The large-scale processed `.h5ad` objects and raw network inference outputs are archived in a permanent Zenodo repository for transparency and reproducibility:

**DOI: 10.5281/zenodo.20715903**
**https://zenodo.org/records/20715903**

## Requirements
This pipeline is built for Python 3.10+ and requires the following core libraries:
`scanpy`, `anndata`, `pyscenic`, `harmonypy`, `liana`, `gseapy`, and `pingouin`.

## License
This project is licensed under the Creative Commons Attribution 4.0 International (CC BY 4.0) license.
