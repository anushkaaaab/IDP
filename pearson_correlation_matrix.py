import pandas as pd
import numpy as np

# Load clean expression subset matrix (samples x genes)
clean_subset = pd.read_csv("clean_expression_subset.csv", index_col=0)

# 1. Compute full Pearson correlation matrix
pearson_matrix = clean_subset.corr(method='pearson')

# Save raw correlation matrix to CSV
pearson_matrix.to_csv("pearson_correlation_matrix.csv")
print("Saved Pearson correlation matrix of shape:", pearson_matrix.shape)

# 2. Apply threshold to inspect edge count
THRESHOLD = 0.65

# Create a explicit WRITABLE copy of the array values
adj_values = pearson_matrix.to_numpy(copy=True)

# Zero out diagonal (self-correlations)
np.fill_diagonal(adj_values, 0)

# Apply thresholding
adj_values[np.abs(adj_values) < THRESHOLD] = 0

# Convert back to DataFrame
adj_matrix = pd.DataFrame(adj_values, index=pearson_matrix.index, columns=pearson_matrix.columns)

# Calculate retained edges
num_edges = np.count_nonzero(adj_values) // 2
total_possible = pearson_matrix.shape[0] * (pearson_matrix.shape[0] - 1) // 2

print(f"Total possible gene pairs: {total_possible}")
print(f"Significant co-expression edges retained at |r| >= {THRESHOLD}: {num_edges}")