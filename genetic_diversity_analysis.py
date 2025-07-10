"""
Complete genetic diversity analysis script.

Before running:
- Make sure you have installed necessary packages:
    pip install numpy pandas scipy matplotlib seaborn scikit-allel

Adjust:
- 'vcf_path' to point to your VCF or genotype data.
- Or replace the example genotype loading with your actual data loading.

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram
import allel  # scikit-allel for VCF/genotype handling

# === Load genotype data ===
# Example: Load a VCF file (replace with your actual VCF file path)
vcf_path = 'targets.dp20gq50.filtered.vcf'  # Replace with your path

print("Loading VCF and extracting genotypes...")
callset = allel.read_vcf(vcf_path, fields=['samples', 'calldata/GT'])
samples = callset['samples']
genotypes = allel.GenotypeArray(callset['calldata/GT'])

# Convert genotypes to allele counts (0,1,2 for diploid)
geno_ac = genotypes.to_n_alt()  # shape (variants, samples)

# Transpose to (samples, variants) for pairwise calculations
geno_ac_T = geno_ac.T

print(f"Number of samples: {len(samples)}")
print(f"Number of variants: {geno_ac.shape[0]}")

# === Compute pairwise genetic distance matrix ===
# Here we use Euclidean distance between genotype vectors (simple approach)
print("Computing pairwise genetic distances...")

from scipy.spatial.distance import pdist, squareform

# pdist expects 2D array shape=(n_samples, n_features)
distance_vector = pdist(geno_ac_T, metric='euclidean')
distance_matrix = squareform(distance_vector)

# Convert to DataFrame for easy handling
df_distance = pd.DataFrame(distance_matrix, index=samples, columns=samples)

# Save distance matrix
df_distance.to_csv('genetic_distance_matrix.csv')
print("Genetic distance matrix saved as 'genetic_distance_matrix.csv'")

# === Plot heatmap ===
plt.figure(figsize=(14, 12))
sns.heatmap(df_distance, cmap='viridis')
plt.title('Genetic Distance Matrix Heatmap')
plt.tight_layout()
plt.savefig('genetic_distance_heatmap.png')
print("Heatmap saved as 'genetic_distance_heatmap.png'")
plt.show()

# === Hierarchical clustering and dendrogram ===
linkage_matrix = linkage(distance_vector, method='average')

plt.figure(figsize=(15, 7))
dendrogram(linkage_matrix, labels=samples, leaf_rotation=90)
plt.title('Hierarchical Clustering Dendrogram')
plt.tight_layout()
plt.savefig('genetic_distance_dendrogram.png')
print("Dendrogram saved as 'genetic_distance_dendrogram.png'")
plt.show()
