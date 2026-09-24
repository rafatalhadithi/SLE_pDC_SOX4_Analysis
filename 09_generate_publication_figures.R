# Large-Cohort Validation of a Pathogenic SOX4-CXCR4 Network in SLE pDCs

**Repository Status:** Private (Peer-Review Phase)

This repository contains the consolidated 9-step computational pipeline used to process, integrate, and perform network inference on plasmacytoid dendritic cell (pDC) scRNA-seq data from patients with Systemic Lupus Erythematosus (SLE) for the manuscript: *"Large-Cohort Validation of a Pathogenic SOX4-CXCR4 Network in SLE pDCs."*

## Overview
This workflow utilizes `Scanpy`, `pySCENIC`, `LIANA`, and `Harmony` to identify the **SOX4-CXCR4** regulatory axis. The analysis isolates the pathogenic network and leverages large-scale cohorts to mathematically decouple systemic inflammation from tissue-specific homing mechanisms in Lupus Nephritis.

---

## 🚀 Quick Start Guide for Reviewers

To ensure seamless reproducibility across different operating systems, this pipeline uses dynamic relative paths. **Ensure your terminal or IDE working directory is set to the repository root before executing the scripts.**

**1. Clone the Repository**
`git clone https://github.com/YOUR_USERNAME/SLE_pDC_SOX4_Analysis.git`
`cd SLE_pDC_SOX4_Analysis`

**2. Set Up the Data Environment**
The scripts will automatically generate `outputs/` and `figures/` directories when run. Populate the `data/` folder with the required datasets before executing the pipeline:
* **Processed Matrices:** Download `SLE_pDC_Cleaned_For_pySCENIC.h5ad`, `SLE_pDC_pySCENIC_Complete.h5ad`, and `SLE_pDC_Integrated.h5ad` from our Zenodo repository (link below).
* **Global Atlas:** Download the 1.2M cell adult SLE PBMC cohort (`4118e166-34f5-4c1f-9eed-c64b90a3dace.h5ad`) from the CZ CELLxGENE portal.
* **pySCENIC Databases:** Download `hs_hgnc_tfs.txt`, the hg38 motif `.feather` database, and the motif annotation `.tbl` from the cisTarget repository.
* **Validation Cohorts:** Download `GSE135779` and `GSE303481` raw data from NCBI GEO and place them in `data/GSE135779_RAW_Data/` and `data/GSE303481_RAW_Data/` respectively.

**3. Execute the Pipeline**
Run the scripts sequentially from `01` to `09`.

---

## Repository Contents & Exact File I/O Mapping

### `01_preprocess_global_ifn.py`
Initial quality control, pDC isolation, and global Type I IFN volumetric profiling.
* **Input File:** `data/4118e166-34f5-4c1f-9eed-c64b90a3dace.h5ad`
* **Output Artifacts:** `data/SLE_pDC_Cleaned_For_pySCENIC.h5ad`, `outputs/01_Global_Type1_IFN_Production_Summary.csv`, `outputs/SLE_pDC_Top77_Upregulated.csv`, `figures/01_Global_Type1_IFN_Production.pdf`

### `02_pyscenic_grn_inference.py`
Executes machine learning gene regulatory network inference via GRNBoost2 and RcisTarget.
* **Input Files:** `data/SLE_pDC_Cleaned_For_pySCENIC.h5ad`, `data/hs_hgnc_tfs.txt`, `data/hg38__refseq-r80...feather`, `data/motifs-v9-nr...tbl`
* **Output Artifacts:** `data/SLE_pDC_pySCENIC_Complete.h5ad`, `data/IRF7_Target_Genes.txt`, `outputs/Table_S1_Regulons_Raw.csv`

### `03_go_pathway_enrichment.py`
Performs Gene Ontology pathway enrichment testing for the IRF7 regulon.
* **Input File:** `data/IRF7_Target_Genes.txt`
* **Output Artifacts:** `outputs/Table_S2_GO_Enrichment.csv`, `figures/Figure_2_IRF7_GO_Enrichment.pdf`

### `04_epigenetic_atac_profiling.py`
Orthogonal validation querying ENCODE ATAC-seq databases to map physical SOX4 DNA binding motifs at the CXCR4 promoter.
* **Input File:** ENCODE PBMC ATAC-seq data (Auto-fetches via API to `data/Human_PBMC_ATACseq_narrowPeak.bed.gz`)
* **Output Artifacts:** `outputs/04_SOX4_Motifs.csv`, `figures/Figure_S4_ATAC_Motif_Footprint.pdf`

### `05_liana_cell_communication.py`
Maps ligand-receptor interaction networks signaling into the pDC target cells.
* **Input File:** `data/4118e166-34f5-4c1f-9eed-c64b90a3dace.h5ad`
* **Output Artifacts:** `outputs/LIANA_Full_Results.csv`, `outputs/pDC_CXCR4_Interactions.csv`

### `06_statistical_power_and_confounders.py`
Executes donor-level pseudobulk partial correlations (regressing out the IFN signature) and a 1,000-iteration bootstrapping power simulation.
* **Input File:** `data/SLE_pDC_pySCENIC_Complete.h5ad`
* **Output Artifacts:** Console statistical readouts, `figures/Figure_4_Bootstrapping_Power_Analysis.pdf`

### `07_harmony_batch_integration.py`
Executes native batch-effect correction across 242 distinct donor profiles.
* **Input File:** `data/SLE_pDC_pySCENIC_Complete.h5ad`
* **Output Artifacts:** `data/SLE_pDC_Integrated.h5ad`

### `08_cross_cohort_meta_analysis.py`
Independently evaluates the SOX4-CXCR4 axis in an external pediatric cohort (GSE135779) and adult Lupus Nephritis kidney cohort (GSE303481), synthesizing p-values via Fisher’s Combined Probability Test.
* **Input Files:** `data/GSE135779_genes.tsv.gz`, `data/GSE135779_RAW_Data/`, `data/GSE303481_RAW_Data/`
* **Output Artifacts:** Console statistical readouts

### `09_generate_publication_figures.R`
Generates high-resolution publication-grade visual figures utilizing `ggplot2` and `reticulate`.
* **Input Files:** `data/SLE_pDC_Integrated.h5ad`
* **Output Artifacts:** `figures/Supplementary_FigS1_Harmony_Integration.pdf`, `figures/Fig1A_Network_vs_Gene_Validation.pdf`, `figures/Fig1B_Top_SLE_Regulons_Dotplot.pdf`

---

## Data Availability
The large-scale processed `.h5ad` objects and raw network inference outputs are archived in a permanent Zenodo repository for transparency and reproducibility:

**DOI: 10.5281/zenodo.20715903**  
**Confidential Reviewer Access Link (Bypasses Embargo):** [Click Here for Data Access](https://zenodo.org/records/20715903?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjBmMzQ1YmMyLTQyNDYtNDcwOC1iMzM3LWM0YWYwYTJiZWYzNyIsImRhdGEiOnt9LCJyYW5kb20iOiJhYTNhZmJlM2I1YTM5N2ZkOTNmYmM2YTcxMmQ2OGJiNiJ9.7sOJI_Hdk_xrtas08gGL_YkhmVi69BBP4aZpOPv0w3kfutyMJvpBTKOhcRWvGtQoBtUl2BZSk4I6ngrRjBrcrw)

---

## Requirements & Installation

This pipeline is built for Python 3.10+ and R 4.x. 

### 1. Install Python Dependencies
```bash
pip install scanpy pandas numpy scipy matplotlib seaborn pyscenic dask distributed arboreto harmonypy liana gseapy pingouin
```

### 2. Install R Dependencies (For Figure Generation)
Script `09_generate_publication_figures.R` requires R and the following packages. Run this in your R console:

```R
install.packages(c("reticulate", "ggplot2", "dplyr", "patchwork", "viridis"))
```

# License
This project is licensed under the Creative Commons Attribution 4.0 International (CC BY 4.0) license.
