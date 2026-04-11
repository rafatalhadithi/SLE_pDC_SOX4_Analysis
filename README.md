# SLE pDC SOX4-CXCR4 Analysis

**Repository Status:** Private (Peer-Review Phase)

This repository contains the custom Python scripts used for the single-cell RNA-sequencing analysis in our manuscript: *"In Silico Discovery of a Putative SOX4-CXCR4 Network and IRF7-Associated Metabolic Reprogramming in SLE pDCs."*

## Description of Files

Below is a breakdown of the specific scripts used to validate our findings:

* **`validate_kidney.py`**
  * **Purpose:** External validation of the SOX4-CXCR4 correlation within a target organ.
  * **Description:** This script loads raw scRNA-seq matrices from adult Lupus Nephritis kidney biopsies (GSE135779), isolates the infiltrating pDCs, and calculates the Spearman correlation between *SOX4* and *CXCR4* expression to confirm that the network remains active outside of peripheral blood.

* **`reviewer_stats.py`**
  * **Purpose:** Advanced statistical validation for peer review.
  * **Description:** Calculates both donor-level pseudobulk correlations and partial Spearman correlations (regressing out the prominent Type I IFN signature) to prove the SOX4-CXCR4 relationship is robust and independent of general cellular hyperactivation.

* **`run_liana_safe.py`**
  * **Purpose:** Cell-Cell Communication Analysis.
  * **Description:** Utilizes the LIANA framework with consensus resources to map ligand-receptor interactions across the PBMC cohort. This script specifically screens for cells sending signals to pDC CXCR4 receptors, confirming the absence of peripheral targeting.

*(Note: If you uploaded any other scripts, like `run_network.py`, you can just type a quick bullet point for them here!)*

## Data Availability
The massive, heavily processed AnnData (`.h5ad`) matrices, pySCENIC database files, and high-resolution figures utilized by these scripts are too large to host on GitHub. They are safely archived locally and are available from the corresponding author upon reasonable request.
