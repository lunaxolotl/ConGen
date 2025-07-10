import allel
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
from sklearn.decomposition import PCA
import os
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
import matplotlib.colors as mcolors


# Create subfolder "bottleneck" if it doesn't exist
output_dir = "bottleneck"
os.makedirs(output_dir, exist_ok=True)

# === 1. Load VCF ===
vcf_path = "targets.dp20gq50.filtered.vcf"  # Replace with your actual VCF file
callset = allel.read_vcf(vcf_path, fields=["samples", "calldata/GT", "variants/POS"])
gt = allel.GenotypeArray(callset["calldata/GT"])
samples = callset["samples"]
positions = callset["variants/POS"]

# === 2. Load metadata ===
meta = pd.read_csv("AlpineIbex_sampleinfo.csv")  # Must contain 'vcf_sample_id' and 'population'
meta = meta.set_index("vcf_sample_id")
pop_map = meta["population"].to_dict()

# === 3. Map populations to sample indices ===
pop_indices = {}
for pop in meta["population"].unique():
    pop_indices[pop] = [i for i, s in enumerate(samples) if s in meta.index and pop_map[s] == pop]

# === 4. Calculate per-population statistics ===
results = []
for pop, indices in pop_indices.items():
    gt_sub = gt[:, indices]
    if gt_sub.n_variants == 0 or gt_sub.n_samples == 0:
        continue

    ac = gt_sub.count_alleles()
    ho = gt_sub.count_het(axis=1).mean()  # observed heterozygosity

    af = ac.to_frequencies()
    he_per_variant = 1 - np.sum(af**2, axis=1)
    he = np.nanmean(he_per_variant)  # expected heterozygosity

    # Nucleotide diversity (pi)
    sorted_indices = np.argsort(positions)
    positions_sorted = positions[sorted_indices]
    ac_sorted = ac[sorted_indices]
    pi = allel.sequence_diversity(positions_sorted, ac_sorted)

    # Tajima's D - only biallelic variants, using derived allele counts (alternate allele)
    try:
        is_biallelic = ac.is_biallelic()
        ac_bi = ac.compress(is_biallelic, axis=0)
        tajd = allel.tajima_d(ac_bi)
        tajd_mean = np.nanmean(tajd)
    except Exception as e:
        print(f"Tajima's D error in population {pop}: {e}")
        tajd_mean = np.nan

    results.append({
        "population": pop,
        "heterozygosity_obs": ho,
        "heterozygosity_exp": he,
        "nucleotide_diversity": pi,
        "tajimas_D": tajd_mean
    })

df_results = pd.DataFrame(results)
print("\n=== Population Summary Statistics ===")
print(df_results)

# Export population statistics table
df_results.to_csv(os.path.join(output_dir, "population_summary_statistics.csv"), index=False)

# === 5. Load and process genetic distance matrix ===
dist_df = pd.read_csv("genetic_distance_matrix.csv", index_col=0)

def pop_distance(pop1, pop2):
    s1 = meta[meta["population"] == pop1].index
    s2 = meta[meta["population"] == pop2].index
    common = dist_df.index.intersection(s1)
    common_cols = dist_df.columns.intersection(s2)
    if len(common) == 0 or len(common_cols) == 0:
        return np.nan
    return dist_df.loc[common, common_cols].values.mean()

pops = df_results["population"].tolist()
intensity = pd.DataFrame(index=pops, columns=pops, dtype=float)

for i in pops:
    for j in pops:
        intensity.loc[i, j] = pop_distance(i, j)

plt.figure(figsize=(10, 8))
sns.heatmap(intensity, annot=True, fmt=".2f", cmap="viridis", linewidths=0.5)
plt.title("Mean Genetic Distance Between Populations")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "mean_genetic_distance_heatmap.png"))
plt.show()

# Export genetic distance matrix
intensity.to_csv(os.path.join(output_dir, "mean_genetic_distance_matrix.csv"))

# === 6. Calculate pairwise FST between populations ===
fst_results = []
for pop1, pop2 in combinations(pops, 2):
    idx1 = pop_indices[pop1]
    idx2 = pop_indices[pop2]
    if len(idx1) == 0 or len(idx2) == 0:
        continue

    # Combine genotype array for both populations
    combined_idx = np.array(idx1 + idx2)
    gt_sub = gt[:, combined_idx]

    # Calculate allele counts for each subpopulation in combined data
    ac1 = gt_sub.count_alleles(subpop=range(len(idx1)))
    ac2 = gt_sub.count_alleles(subpop=range(len(idx1), len(combined_idx)))

    # Filter biallelic variants in both subpops
    is_bi = ac1.is_biallelic() & ac2.is_biallelic()
    gt_bi = gt_sub.compress(is_bi, axis=0)

    # subpops indices relative to gt_bi samples
    subpops = [np.arange(len(idx1)), np.arange(len(idx1), len(combined_idx))]

    try:
        fst = allel.weir_cockerham_fst(gt_bi, subpops=subpops)
        fst_mean = np.nanmean(fst)
    except Exception as e:
        print(f"FST calculation error for {pop1}-{pop2}: {e}")
        fst_mean = np.nan

    fst_results.append({"pop1": pop1, "pop2": pop2, "FST": fst_mean})

df_fst = pd.DataFrame(fst_results)
print("\n=== Pairwise FST between populations ===")
print(df_fst)

# Export pairwise FST results
df_fst.to_csv(os.path.join(output_dir, "pairwise_fst.csv"), index=False)

# === 7. PCA on genotype data ===
# Get alternate allele counts per variant per sample: shape (variants, samples)
ac = gt.to_n_alt()
ac_T = ac.T  # samples x variants

# Remove variants with zero variance
variant_var = np.var(ac_T, axis=0)
ac_T_filtered = ac_T[:, variant_var > 0]

pca = PCA(n_components=2)
pcs = pca.fit_transform(ac_T_filtered)

plt.figure(figsize=(8, 6))
for pop in pops:
    idx = [i for i, s in enumerate(samples) if pop_map.get(s) == pop]
    plt.scatter(pcs[idx, 0], pcs[idx, 1], label=pop, alpha=0.7)
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.2f}%)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.2f}%)")
plt.title("PCA of Genotype Data")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "pca_genotype_data.png"))
plt.show()

# === 8. Plot He vs Ho ===
plt.figure(figsize=(6, 6))
plt.scatter(df_results["heterozygosity_exp"], df_results["heterozygosity_obs"], c='blue')
plt.plot([0, 1], [0, 1], 'k--')
plt.xlabel("Expected Heterozygosity (He)")
plt.ylabel("Observed Heterozygosity (Ho)")
plt.title("He vs Ho per Population")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "he_vs_ho.png"))
plt.show()

# === 9. Calculate pairwise relatedness matrix (IBS) and plot heatmap ===
print("\n=== Calculating pairwise relatedness matrix (IBS) ===")

# Get alternate allele counts (variants x samples)
ac = gt.to_n_alt()
ac_T = ac.T  # samples x variants

n_samples = ac_T.shape[0]
ibs_matrix = np.zeros((n_samples, n_samples))

for i in range(n_samples):
    for j in range(i, n_samples):
        valid = ~(np.isnan(ac_T[i]) | np.isnan(ac_T[j]))
        n_valid = valid.sum()
        if n_valid == 0:
            ibs = np.nan
        else:
            # IBS similarity scaled between 0 and 1
            ibs = (2 - np.abs(ac_T[i, valid] - ac_T[j, valid])).sum() / (2 * n_valid)
        ibs_matrix[i, j] = ibs
        ibs_matrix[j, i] = ibs

ibs_df = pd.DataFrame(ibs_matrix, index=samples, columns=samples)
ibs_df.to_csv(os.path.join(output_dir, "pairwise_relatedness_matrix.csv"))

plt.figure(figsize=(12, 10))
sns.heatmap(ibs_df, cmap="coolwarm", center=0.5)
plt.title("Pairwise Relatedness (IBS) Matrix")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "pairwise_relatedness_heatmap.png"))
plt.show()

# === 10. Hierarchical clustering dendrogram with population color labels ===

# Step 1: Fill NaNs in intensity with max+1 (large distance)
max_val = intensity.max().max()
intensity_filled = intensity.fillna(max_val + 1)

# Step 2: Symmetrize the matrix by averaging with its transpose
intensity_sym = (intensity_filled + intensity_filled.T) / 2

# Step 3: Set diagonal to zero (distance to self)
np.fill_diagonal(intensity_sym.values, 0)

# Step 4: Convert symmetrized matrix to condensed distance matrix for clustering
condensed_dist = squareform(intensity_sym.values)

# Step 5: Perform hierarchical clustering (average linkage)
Z = linkage(condensed_dist, method='average')

# Step 6: Create a color palette and dict for populations
pop_colors = sns.color_palette("tab10", len(pops))
color_dict = dict(zip(pops, pop_colors))

# Step 7: Plot dendrogram with colored population labels
plt.figure(figsize=(12, 6))
dendro = dendrogram(Z, labels=pops, leaf_rotation=90)

ax = plt.gca()
xlbls = ax.get_xmajorticklabels()
for lbl in xlbls:
    pop = lbl.get_text()
    lbl.set_color(mcolors.to_hex(color_dict[pop]))

plt.title("Hierarchical Clustering of Populations Based on Genetic Distance")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "hierarchical_clustering_dendrogram.png"))
plt.show()

