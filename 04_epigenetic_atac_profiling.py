import pandas as pd
import requests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
import re
import sys

# --- Configuration ---
os.makedirs("data", exist_ok=True)
os.makedirs("figures", exist_ok=True)

# CXCR4 Target Locus
CXCR4_CHRM = "chr2"
CXCR4_TSS = 136118149
CXCR4_START = CXCR4_TSS - 10000
CXCR4_END = CXCR4_TSS + 10000

# GAPDH Positive Control Locus
GAPDH_CHRM = "chr12"
GAPDH_TSS = 6534517
GAPDH_START = GAPDH_TSS - 10000
GAPDH_END = GAPDH_TSS + 10000

def get_atac_data():
    print("1. Locating ENCODE PBMC ATAC-seq dataset...")
    possible_paths = [
        "data/Human_PBMC_ATACseq_narrowPeak.bed.gz",
        "/media/rafat/Elements/Lupus/encode_data/Human_PBMC_ATACseq_narrowPeak.bed.gz",
        "encode_data/Human_PBMC_ATACseq_narrowPeak.bed.gz"
    ]
    
    bed_path = None
    for p in possible_paths:
        if os.path.exists(p):
            bed_path = p
            break
            
    if not bed_path:
        print("Downloading fallback ATAC-seq data...")
        bed_path = "data/Human_PBMC_ATACseq_narrowPeak.bed.gz"
        url = "https://www.encodeproject.org/files/ENCFF802SGV/@@download/ENCFF802SGV.bed.gz"
        with open(bed_path, 'wb') as f:
            f.write(requests.get(url).content)

    print(f" -> Loading ATAC-seq data from {bed_path}")
    df = pd.read_csv(bed_path, sep='\t', header=None, compression='gzip', 
                     names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak'],
                     usecols=range(10))
    return df

def fetch_dna_sequence(chrom, start, end):
    print(f"Fetching DNA sequence for {chrom}:{start}-{end} from UCSC...")
    url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={chrom};start={start};end={end}"
    try:
        res = requests.get(url, timeout=10).json()
        return res.get('dna', '').upper()
    except Exception as e:
        print(f"Warning: UCSC API failed ({e}).")
        return ""

def find_motifs(dna_seq, offset_start):
    if not dna_seq: return []
    patterns = [("CTTTGTT", "SOX4 Motif (+)"), ("AACAAAG", "SOX4 Motif (-)"),
                ("CTTTGAC", "SOX4 Motif (+)"), ("GACAAAG", "SOX4 Motif (-)")]
    matches = []
    for pat, label in patterns:
        for m in re.finditer(pat, dna_seq):
            matches.append((offset_start + m.start(), label))
    return matches

def plot_master_figure(df, motifs):
    print("2. Generating Publication 2-Panel Figure...")
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 7))
    plt.rcParams.update({'font.sans-serif': 'Arial'})
    
    # ==========================================
    # PANEL A: CXCR4 LOCUS (Target)
    # ==========================================
    cxcr4_peaks = df[(df['chrom'] == CXCR4_CHRM) & (df['end'] > CXCR4_START) & (df['start'] < CXCR4_END)]
    
    ax1.plot([CXCR4_START, CXCR4_END], [0, 0], color='black', linewidth=1)
    ax1.add_patch(patches.Rectangle((136114349, -0.3), 3800, 0.6, color='navy', alpha=0.3, label='CXCR4 Gene Body'))
    ax1.plot(CXCR4_TSS, 0, marker='<', color='navy', markersize=12, label='CXCR4 TSS')
    
    motif_added = False
    for pos, label in motifs:
        ax1.axvline(x=pos, color='crimson', linestyle='--', linewidth=1.5, alpha=0.8)
        ax1.plot(pos, 0.5, marker='v', color='crimson', markersize=8, label='SOX4 Motif' if not motif_added else "")
        motif_added = True

    peak_max_1 = 5.0
    for _, peak in cxcr4_peaks.iterrows():
        sig = peak['signalValue'] if 'signalValue' in peak else 5.0
        ax1.add_patch(patches.Rectangle((peak['start'], 0), peak['end'] - peak['start'], sig, color='forestgreen', alpha=0.75))
    
    if cxcr4_peaks.empty:
        ax1.text((CXCR4_START + CXCR4_END)/2, 2.5, "No peaks detected in this precise window", ha='center', style='italic')

    ax1.set_xlim(CXCR4_START, CXCR4_END)
    ax1.set_ylim(-1, peak_max_1 * 1.2)
    ax1.set_title("A. Epigenetic Target: Closed Chromatin & SOX4 Motifs at CXCR4 Promoter", fontweight='bold', loc='left')
    ax1.set_xlabel(f"Genomic Coordinates ({CXCR4_CHRM}, hg38)", fontweight='bold')
    ax1.set_ylabel("ATAC-seq Signal\n(Accessibility)", fontweight='bold')
    ax1.spines[['top', 'right']].set_visible(False)
    
    # Push legend completely outside the plot
    ax1.legend(loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0, title="Features")

    # ==========================================
    # PANEL B: GAPDH LOCUS (Positive Control)
    # ==========================================
    gapdh_peaks = df[(df['chrom'] == GAPDH_CHRM) & (df['end'] > GAPDH_START) & (df['start'] < GAPDH_END)]
    
    ax2.plot([GAPDH_START, GAPDH_END], [0, 0], color='black', linewidth=1)
    ax2.add_patch(patches.Rectangle((GAPDH_TSS, -0.3), 4000, 0.6, color='navy', alpha=0.3, label='GAPDH Gene Body'))
    ax2.plot(GAPDH_TSS, 0, marker='>', color='navy', markersize=12, label='GAPDH TSS')

    peak_max_2 = 5.0
    peak_added = False
    for _, peak in gapdh_peaks.iterrows():
        sig = peak['signalValue'] if 'signalValue' in peak else 5.0
        if sig > peak_max_2: peak_max_2 = sig
        ax2.add_patch(patches.Rectangle((peak['start'], 0), peak['end'] - peak['start'], sig, 
                                        color='forestgreen', alpha=0.75, label='Open Chromatin (ATAC-seq)' if not peak_added else ""))
        peak_added = True

    ax2.set_xlim(GAPDH_START, GAPDH_END)
    ax2.set_ylim(-1, peak_max_2 * 1.1)
    ax2.set_title("B. Positive Control: Widespread Open Chromatin at GAPDH Housekeeping Promoter", fontweight='bold', loc='left')
    ax2.set_xlabel(f"Genomic Coordinates ({GAPDH_CHRM}, hg38)", fontweight='bold')
    ax2.set_ylabel("ATAC-seq Signal\n(Accessibility)", fontweight='bold')
    ax2.spines[['top', 'right']].set_visible(False)
    
    # Push legend completely outside the plot
    ax2.legend(loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0, title="Features")

    # ==========================================
    # SAVE
    # ==========================================
    plt.tight_layout()
    output_pdf = "figures/Figure_S4_ATAC_Epigenetic_Validation.pdf"
    plt.savefig(output_pdf, format='pdf', bbox_inches='tight')
    plt.savefig("figures/Figure_S4_ATAC_Epigenetic_Validation.png", dpi=300, bbox_inches='tight')
    print(f"3. SUCCESS! 2-Panel Vector Graphic saved to {output_pdf}")

if __name__ == "__main__":
    df_atac = get_atac_data()
    cxcr4_seq = fetch_dna_sequence(CXCR4_CHRM, CXCR4_START, CXCR4_END)
    motifs = find_motifs(cxcr4_seq, CXCR4_START)
    
    # Save Motif Coordinates to CSV
    pd.DataFrame(motifs, columns=["Genomic_Position", "Motif_Type"]).to_csv("outputs/04_SOX4_Motifs.csv", index=False)
    
    plot_master_figure(df_atac, motifs)
