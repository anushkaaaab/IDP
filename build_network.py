import pandas as pd
import numpy as np
import os

def generate_sample_data(filename="clean_expression_data.csv", n_genes=50, n_samples=100):
    print(f"Generating synthetic expression data in '{filename}'...")

    np.random.seed(42)

    genes = [f"GENE_{i}" for i in range(1, n_genes - 4)] + [
        "TP53", "ESR1", "ERBB2", "BRCA1", "BRCA2"
    ]

    data = np.random.randn(len(genes), n_samples)

    for i in range(15):
        data[genes.index("TP53")] += data[i] * 0.8
        data[genes.index("ESR1")] += data[i + 10] * 0.75
        data[genes.index("ERBB2")] += data[i + 20] * 0.7

    df = pd.DataFrame(
        data,
        index=genes,
        columns=[f"Sample_{i}" for i in range(1, n_samples + 1)]
    )

    df.to_csv(filename, index_label="gene")
    print(f"Created expression data: {df.shape[0]} genes, {df.shape[1]} samples")


def build_network_edges(
    expression_csv="clean_expression_data.csv",
    output_csv="network_edges.csv",
    threshold=0.3
):
    if not os.path.exists(expression_csv) or os.path.getsize(expression_csv) == 0:
        generate_sample_data(expression_csv)

    try:
        df = pd.read_csv(expression_csv, index_col=0)
    except pd.errors.EmptyDataError:
        print("Input CSV is empty. Regenerating sample data.")
        generate_sample_data(expression_csv)
        df = pd.read_csv(expression_csv, index_col=0)

    if df.empty or df.shape[1] == 0:
        print("No expression data found.")
        return

    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(axis=0, how="all").dropna(axis=1, how="all")

    print(f"Loaded expression data: {df.shape[0]} genes, {df.shape[1]} samples")

    if df.shape[0] < 2 or df.shape[1] < 2:
        print("At least 2 genes and 2 samples are required.")
        return

    corr_matrix = df.T.corr(method="pearson")

    edges = []

    genes = corr_matrix.columns

    for i in range(len(genes)):
        for j in range(i + 1, len(genes)):
            weight = corr_matrix.iloc[i, j]

            if pd.notna(weight) and abs(weight) >= threshold:
                edges.append({
                    "gene_1": str(genes[i]),
                    "gene_2": str(genes[j]),
                    "weight": float(weight)
                })

    edge_df = pd.DataFrame(edges, columns=["gene_1", "gene_2", "weight"])
    edge_df.to_csv(output_csv, index=False)

    print(f"Highest absolute correlation: {corr_matrix.abs().where(~np.eye(len(genes), dtype=bool)).max().max():.4f}")
    print(f"Threshold: {threshold}")
    print(f"Created '{output_csv}' with {len(edge_df)} network edges.")

    if edge_df.empty:
        print("No gene pairs passed the correlation threshold.")
        print("Try threshold=0.3 to test whether the threshold is too strict.")
    else:
        print(edge_df.head(10))


if __name__ == "__main__":
    build_network_edges()