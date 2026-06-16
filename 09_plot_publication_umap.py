import scanpy as sc
import matplotlib.pyplot as plt

print("1. Loading integrated dataset...")
adata = sc.read_h5ad("SLE_pDC_Integrated.h5ad")

print("2. Setting up publication-grade figure layout...")
# Create a wide figure with 3 side-by-side subplots
fig, axs = plt.subplots(1, 3, figsize=(18, 5))

print("3. Generating plots...")
# Panel 1: donor_id (Turning OFF the massive legend so it fits cleanly)
sc.pl.umap(adata, color='donor_id', ax=axs[0], legend_loc=None, show=False, 
           title='Harmony Integration\n(242 Donors Blended)')

# Panel 2: disease state (Keeping standard legend)
sc.pl.umap(adata, color='disease', ax=axs[1], show=False, 
           title='Disease State')

# Panel 3: SOX4 expression (Keeping colorbar, using viridis for high contrast)
sc.pl.umap(adata, color='SOX4', ax=axs[2], color_map='viridis', show=False, 
           title='SOX4 Expression')

print("4. Formatting and saving to PDF...")
# Ensure spacing is perfect so titles and legends don't overlap
plt.tight_layout()

# Save as a high-resolution vector PDF for the manuscript
output_file = "Figure_Harmony_UMAP_Publication.pdf"
plt.savefig(output_file, dpi=300, bbox_inches='tight')

print(f"\n>>> Masterpiece created! Saved as: {output_file} <<<")
