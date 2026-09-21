#!/usr/bin/env python3
"""
Script: 04_epigenetic_atac_profiling.py
Purpose: Queries ENCODE ATAC-seq data to map physical SOX4 footprints at the CXCR4 promoter.
"""
import os
import requests
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import re

BASE_DIR = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
DATA_DIR = os.path.join(BASE_DIR, "data")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")
FIGURES_DIR = os.path.join(BASE_DIR, "figures")

os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

CHRM = "chr2"
CXCR4_TSS = 136118149 
WINDOW = 10000 
START, END = CXCR4_TSS - WINDOW, CXCR4_TSS + WINDOW

def fetch_pbmc_atac():
    print("Fetching ENCODE Primary Human PBMC ATAC-seq dataset (ENCFF802SGV)...")
    bed_url = "https://www.encodeproject.org/files/ENCFF802SGV/@@download/ENCFF802SGV.bed.gz"
    local_bed = os.path.join(DATA_DIR, "Human_PBMC_ATACseq_narrowPeak.bed.gz")
    
    if not os.path.exists(local_bed):
        res = requests.get(bed_url)
        with open(local_bed, 'wb') as f:
            f.write(res.content)
            
    df = pd.read_csv(local_bed, sep='\t', header=None, compression='gzip', usecols=range(10),
                     names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'sig', 'pval', 'qval', 'peak'])
    return df

def fetch_cxcr4_dna():
    print("Fetching hg38 DNA sequence from UCSC API...")
    url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={CHRM};start={START};end={END}"
    return requests.get(url).json().get('dna', '').upper()

if __name__ == "__main__":
    atac_df = fetch_pbmc_atac()
    dna_seq = fetch_cxcr4_dna()
    
    # Scan for SOX4 motifs
    motifs = []
    for pattern, label in [("CTTTGTT", "SOX4 (+)"), ("AACAAAG", "SOX4 (-)"), ("CTTTGAC", "SOX4 (+)"), ("GACAAAG", "SOX4 (-)")]:
        for m in re.finditer(pattern, dna_seq):
            motifs.append((START + m.start(), label))
            
    pd.DataFrame(motifs, columns=["Genomic_Position", "Motif_Type"]).to_csv(os.path.join(OUTPUT_DIR, "04_SOX4_Motifs.csv"), index=False)
    print(f"Mapped {len(motifs)} SOX4 motifs. Generating Genomic Track Plot...")

    locus_peaks = atac_df[(atac_df['chrom'] == CHRM) & (atac_df['end'] > START) & (atac_df['start'] < END)]
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 6), sharex=True, gridspec_kw={'height_ratios': [3, 1.2]})
    
    ax1.plot([START, END], [0, 0], color='black', linewidth=1)
    if locus_peaks.empty:
        ax1.text((START + END) / 2, 2, "No Open Chromatin Peaks in Healthy PBMCs", ha='center', style='italic')
    
    ax1.set_ylabel("ATAC-seq Signal", fontweight='bold')
    ax1.set_title("Epigenetic Validation: Chromatin Accessibility & SOX4 Motifs at CXCR4 (hg38)", fontweight='bold')
    ax1.set_ylim(-0.5, 10)
    ax1.spines[['top', 'right']].set_visible(False)

    ax2.plot([START, END], [0, 0], color='black', linewidth=1)
    ax2.add_patch(patches.Rectangle((136114349, -0.3), 3800, 0.6, color='navy', alpha=0.35))
    ax2.plot(CXCR4_TSS, 0, marker='<', color='navy', markersize=12, label='CXCR4 TSS')
    
    for pos, label in motifs:
        ax2.axvline(x=pos, color='crimson', linestyle='--', alpha=0.8)
        
    ax2.set_xlabel(f"Genomic Coordinates ({CHRM})", fontweight='bold')
    ax2.set_yticks([])
    ax2.spines[['top', 'right']].set_visible(False)
    
    plt.savefig(os.path.join(FIGURES_DIR, "Figure_S4_ATAC_Motif_Footprint.pdf"), format='pdf', bbox_inches='tight')
    print("SUCCESS: Epigenetic track saved.")
