import scanpy as sc
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns

print("1. Loading the Blood SLE pDC dataset...")
adata = sc.read_h5ad("SLE_pDC_Cleaned_For_pySCENIC.h5ad")

sox4_full = adata[:, 'SOX4'].X.toarray().flatten()
cxcr4_full = adata[:, 'CXCR4'].X.toarray().flatten()

print("2. Running Bootstrapping Simulation (1,000 iterations at n=127)...")
n_iterations = 1000
sample_size = 127

significant_count = 0
r_values = []

for i in range(n_iterations):
    # Randomly select exactly 127 cells from the blood dataset
    random_indices = np.random.choice(len(sox4_full), size=sample_size, replace=False)
    
    sox4_sample = sox4_full[random_indices]
    cxcr4_sample = cxcr4_full[random_indices]
    
    # Calculate correlation for this tiny subset
    corr, pval = stats.spearmanr(sox4_sample, cxcr4_sample)
    r_values.append(corr)
    
    # Check if it would still be considered statistically significant
    if pval < 0.05 and corr > 0:
        significant_count += 1

power_percentage = (significant_count / n_iterations) * 100

print("\n========================================================")
print("   BOOTSTRAPPING POWER ANALYSIS RESULTS (n=127)")
print("========================================================")
print(f"Times the correlation remained significant: {significant_count} / 1000")
print(f"Statistical Power to detect the network:    {power_percentage:.1f}%")

# Generate a histogram to show the reviewer
plt.figure(figsize=(8, 5))
sns.histplot(r_values, bins=30, kde=True, color='purple')
plt.axvline(x=0.28, color='red', linestyle='--', label='Original Blood Effect Size (R=0.28)')
plt.title('Simulated SOX4-CXCR4 Correlations in Blood (n=127 cells)', fontsize=14)
plt.xlabel('Spearman R Value', fontsize=12)
plt.ylabel('Frequency', fontsize=12)
plt.legend()
plt.grid(axis='y', linestyle='--', alpha=0.7)

plt.savefig('Bootstrapping_Power_Analysis.png', dpi=300, bbox_inches='tight')
print(">>> Saved simulation plot as 'Bootstrapping_Power_Analysis.png' <<<")
