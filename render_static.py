import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import os

# Import the network generator from build_network.py
from build_network import build_network_edges

def generate_static_plot(csv_path="network_edges.csv", output_img="static_network.png"):
    # If network_edges.csv does not exist or is empty, generate it automatically
    if not os.path.exists(csv_path) or os.path.getsize(csv_path) == 0:
        print(f"'{csv_path}' not found. Generating it now via build_network...")
        build_network_edges(output_csv=csv_path)

    df = pd.read_csv(csv_path)
    G = nx.Graph()
    for _, row in df.iterrows():
        G.add_edge(str(row['gene_1']), str(row['gene_2']), weight=float(row['weight']))

    plt.figure(figsize=(12, 10))
    pos = nx.spring_layout(G, k=0.15, seed=42)

    # Highlight known breast cancer driver genes
    key_drivers = ["TP53", "ESR1", "ERBB2", "BRCA1", "BRCA2", "PGR"]
    node_colors = ['#FF4B4B' if node in key_drivers else '#1C83E1' for node in G.nodes()]
    node_sizes = [G.degree(node) * 60 + 100 for node in G.nodes()]

    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, alpha=0.85)
    nx.draw_networkx_edges(G, pos, alpha=0.2, edge_color='#888888')
    
    # Label top hubs and driver genes
    labels = {node: node for node in G.nodes() if node in key_drivers or G.degree(node) > 5}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=9, font_weight='bold')

    plt.title("Static Baseline Gene Regulatory Network (Review III)", fontsize=14, fontweight='bold')
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output_img, dpi=300)
    plt.close()
    print(f"Static visualization saved to {output_img}")

if __name__ == "__main__":
    generate_static_plot()