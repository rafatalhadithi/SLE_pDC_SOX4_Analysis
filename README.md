# SLE pDC SOX4-CXCR4 Analysis

**Repository Status:** Private (Peer-Review Phase)

This repository contains the custom Python scripts used for the single-cell RNA-sequencing analysis in our manuscript: *"In Silico Discovery of a Putative SOX4-CXCR4 Network and IRF7-Associated Metabolic Reprogramming in SLE pDCs."*

## Description of Files

Below is a breakdown of the specific scripts used to execute our pipeline and validate our findings:

* **`run_network.py`**
  * **Purpose:** Gene Regulatory Network (GRN) Inference.
  * **Description:** This script initializes a controlled Dask cluster to run the GRNBoost2 machine learning algorithm (the first step of the pySCENIC pipeline). It processes the cleaned SLE pDC expression matrix against a database of human transcription factors to infer co-expression modules and regulatory adjacencies.

* **`run_go_enrichment.py`**
  * **Purpose:** Pathway Enrichment Analysis.
  * **Description:** Utilizes the `gseapy` library to connect to the Enrichr database. It takes the downstream target genes of the identified master regulator (*IRF7*) and calculates statistically significant Gene Ontology (GO) Biological Process pathways, generating the resulting bar plots (e.g., highlighting oxidative phosphorylation).

* **`run_liana_safe.py`**
  * **Purpose:** Cell-Cell Communication Analysis.
  * **Description:** Utilizes the LIANA framework with consensus resources to map ligand-receptor interactions across the PBMC cohort. This script specifically screens for cells sending signals to pDC *CXCR4* receptors, assessing the presence or absence of peripheral targeting (e.g., via *CXCL12*).

* **`validate_kidney.py`**
  * **Purpose:** Target Organ Validation (Lupus Nephritis).
  * **Description:** This script loads raw scRNA-seq matrices from 56 adult Lupus Nephritis kidney biopsies (GSE135779). It digitally isolates the tissue-infiltrating pDCs using gold-standard markers (*CLEC4C*, *LILRA4*) and calculates the Spearman correlation between *SOX4* and *CXCR4* expression to confirm the network remains active inside the inflamed organ.

* **`reviewer_stats.py`**
  * **Purpose:** Advanced Statistical Validation for Peer Review.
  * **Description:** Calculates both donor-level pseudobulk correlations and partial Spearman correlations (regressing out the prominent Type I IFN signature score). This proves that the *SOX4-CXCR4* regulatory relationship is statistically robust, reproducible across donors, and independent of general cellular hyperactivation.

## Data Availability
The massive, heavily processed AnnData (`.h5ad`) matrices, pySCENIC database files, and high-resolution figures utilized by these scripts are too large to host on GitHub. They are safely archived locally and are available from the corresponding author upon reasonable request.
