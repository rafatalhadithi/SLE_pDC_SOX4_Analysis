import gseapy as gp
import matplotlib.pyplot as plt

if __name__ == '__main__':
    print("1. Loading IRF7 target genes...")
    with open("IRF7_Target_Genes.txt", "r") as f:
        # Read the genes and remove any blank spaces
        genes = [line.strip() for line in f if line.strip()]
    
    print(f"   Successfully loaded {len(genes)} genes.")

    print("2. Running Gene Ontology (GO) Enrichment Analysis...")
    print("   (Connecting to the Enrichr GO_Biological_Process database...)")
    
    # We query the massive Enrichr database for Human biological processes
    enr = gp.enrichr(gene_list=genes,
                     gene_sets=['GO_Biological_Process_2023'],
                     organism='human',
                     outdir=None)

    print("3. Extracting the Top 10 Biological Pathways...")
    # Sort the results by statistical significance (Adjusted P-value)
    results = enr.results
    top_10 = results.sort_values(by='Adjusted P-value').head(10)
    
    print("\n=========================================")
    print("🔥 TOP 10 BIOLOGICAL PATHWAYS DRIVEN BY IRF7 🔥")
    for i, row in top_10.iterrows():
        # Clean up the pathway name for easy reading
        term = row['Term'].split(' (GO:')[0] 
        pval = row['Adjusted P-value']
        print(f"  - {term} (p-adj: {pval:.2e})")
    print("=========================================\n")

    print("4. Generating Pathway Bar Plot...")
    # Let gseapy save the file directly using the 'ofname' parameter
    ax = gp.barplot(enr.results,
                    column="Adjusted P-value",
                    title="IRF7 Target Genes: GO Biological Process",
                    top_term=10,
                    figsize=(8, 6),
                    color='darkred',
                    ofname="5_IRF7_GO_Enrichment.png") 
    
    print("SUCCESS: Check your folder for '5_IRF7_GO_Enrichment.png'!")
    
    