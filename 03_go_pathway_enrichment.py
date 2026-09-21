#!/usr/bin/env python3
import os, sys, gseapy as gp, matplotlib.pyplot as plt, seaborn as sns, numpy as np

BASE_DIR = os.getcwd()
DATA_DIR, OUTPUT_DIR, FIGURES_DIR = [os.path.join(BASE_DIR, d) for d in ["data", "outputs", "figures"]]

genes_file = os.path.join(DATA_DIR, "IRF7_Target_Genes.txt")

if __name__ == '__main__':
    with open(genes_file, "r") as f:
        genes = [line.strip() for line in f if line.strip()]
        
    print(f"Running GO Enrichment on {len(genes)} IRF7 genes...")
    enr = gp.enrichr(gene_list=genes, gene_sets=['GO_Biological_Process_2023'], organism='human')
    
    results = enr.results.sort_values(by='Adjusted P-value')
    results.to_csv(os.path.join(OUTPUT_DIR, "Table_S2_GO_Enrichment.csv"), index=False)
    
    top_10 = results.head(10).copy()
    top_10['Clean_Term'] = top_10['Term'].apply(lambda x: x.split(' (GO:')[0])
    top_10['-log10(p-adj)'] = -np.log10(top_10['Adjusted P-value'])
    
    plt.figure(figsize=(9, 6))
    plt.rcParams.update({'font.sans-serif': 'Arial'})
    sns.set_style("white")

    sns.barplot(data=top_10, x='-log10(p-adj)', y='Clean_Term', hue='-log10(p-adj)', palette='flare_r', dodge=False, legend=False, width=0.65)
    plt.title("IRF7 Target Genes: GO Biological Process", fontsize=15, fontweight='bold', pad=15)
    plt.xlabel(r"$-\log_{10}$ (Adjusted P-value)", fontsize=12, fontweight='bold')
    plt.ylabel("")
    sns.despine(left=True)
    plt.tight_layout()

    plt.savefig(os.path.join(FIGURES_DIR, "Figure_2_IRF7_GO_Enrichment.pdf"), format='pdf', bbox_inches='tight')
