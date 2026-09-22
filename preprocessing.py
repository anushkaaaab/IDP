import GEOparse
import pandas as pd
import numpy as np
import networkx as nx

print("Downloading / loading GSE3494 via GEOparse...")
gse = GEOparse.get_GEO("GSE3494", destdir="./geo_data")

# Extract the expression matrix (rows: probes, columns: GSM sample IDs)
expression_df = gse.pivot_samples('VALUE')
print(f"Raw expression matrix loaded: {expression_df.shape[0]} probes across {expression_df.shape[1]} samples.")

gpl_key = list(gse.gpls.keys())[0]
platform = gse.gpls[gpl_key]
annotation_df = platform.table[['ID', 'Gene Symbol']].dropna()

# Clean gene symbols
annotation_df['Gene Symbol'] = annotation_df['Gene Symbol'].apply(lambda x: str(x).split('///')[0].strip())
annotation_df = annotation_df[annotation_df['Gene Symbol'] != '']

# Map probes to symbols
expression_df = expression_df.merge(annotation_df, left_index=True, right_on='ID', how='inner')

# When multiple probes map to the exact same gene, resolve duplicates by median expression
expression_by_gene = expression_df.groupby('Gene Symbol').median(numeric_only=True)
if 'ID' in expression_by_gene.columns:
    expression_by_gene = expression_by_gene.drop(columns=['ID'])

# Transpose to obtain standard interface format: (samples x genes)
clean_matrix = expression_by_gene.T
print(f"Mapped and cleaned matrix: {clean_matrix.shape[0]} samples x {clean_matrix.shape[1]} genes.")

if clean_matrix.values.max() > 100:
    clean_matrix = np.log2(clean_matrix + 1)

TOP_N_GENES = 200
KEY_DRIVERS = ['TP53', 'ESR1', 'ERBB2', 'BRCA1', 'BRCA2', 'PGR']

gene_variances = clean_matrix.var(axis=0)
top_genes = gene_variances.nlargest(TOP_N_GENES).index.tolist()

for gene in KEY_DRIVERS:
    if gene in clean_matrix.columns and gene not in top_genes:
        top_genes.append(gene)

subset_matrix = clean_matrix[top_genes]
print(f"Subset expression matrix created: {subset_matrix.shape[1]} genes selected.")

# Save clean matrix for records / Student B
subset_matrix.to_csv("clean_expression_subset.csv")

print("Computing Pearson correlation matrix...")
corr_matrix = subset_matrix.corr(method='pearson')

THRESHOLD = 0.65
adj_matrix = corr_matrix.copy()
np.fill_diagonal(adj_matrix.values, 0)
adj_matrix[adj_matrix.abs() < THRESHOLD] = 0

G = nx.from_pandas_adjacency(adj_matrix.abs())
# Remove isolated nodes (nodes with degree 0 after thresholding)
G.remove_nodes_from(list(nx.isolates(G)))

print(f"Constructed Baseline Graph: {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")

# Hub gene identification using Degree Centrality
degree_cent = nx.degree_centrality(G)
top_hubs = sorted(degree_cent.items(), key=lambda x: x[1], reverse=True)[:10]

print("\n--- Top 10 Hub Genes by Degree Centrality ---")
for gene, score in top_hubs:
    print(f"{gene}: {score:.4f}")

edges_data = []
for u, v, data in G.edges(data=True):
    # Retrieve original signed correlation weight
    weight = corr_matrix.loc[u, v]
    edges_data.append({'gene_1': u, 'gene_2': v, 'weight': round(weight, 4)})

edges_df = pd.DataFrame(edges_data)
edges_df.to_csv("network_edges.csv", index=False)

# Export hub ranking summary
hubs_df = pd.DataFrame(top_hubs, columns=['gene_symbol', 'degree_centrality'])
hubs_df.to_csv("top_hubs.csv", index=False)

print("\nPipeline execution complete!")
print("Outputs saved: 'clean_expression_subset.csv', 'network_edges.csv', and 'top_hubs.csv'.")