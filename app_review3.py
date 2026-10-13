import streamlit as st
import pandas as pd
import networkx as nx
import streamlit.components.v1 as components
from pyvis.network import Network

st.set_page_config(page_title="GRN Perturbation Simulator", layout="wide")

st.title("In-Silico Gene Knockout & Perturbation Simulator")
st.caption("Breast Cancer (GSE3494) Co-Expression Network Analysis — Phase 3 / Review III")

@st.cache_data
def load_graph(csv_path):
    df = pd.read_csv(csv_path)
    G = nx.Graph()
    for _, row in df.iterrows():
        G.add_edge(str(row['gene_1']), str(row['gene_2']), weight=float(row['weight']))
    return G

def render_pyvis_network(graph, highlight_genes=[]):
    net = Network(height="450px", width="100%", bgcolor="#ffffff", font_color="#333333")
    for node in graph.nodes():
        color = "#FF4B4B" if node in highlight_genes else "#1C83E1"
        size = graph.degree(node) * 3 + 10
        net.add_node(node, label=node, color=color, size=size)
    for u, v in graph.edges():
        net.add_edge(u, v)
    net.toggle_physics(True)
    return net.generate_html()

try:
    base_graph = load_graph("network_edges.csv")
except Exception:
    st.error("Missing 'network_edges.csv'. Please run 'python build_network.py' first.")
    st.stop()

# Sidebar Setup
st.sidebar.header("Perturbation Controls")
all_genes = sorted(list(base_graph.nodes()))
default_seed = "TP53" if "TP53" in all_genes else (all_genes[0] if all_genes else None)

target_genes = st.sidebar.multiselect(
    "Select Target Gene(s) for Structural Knockout:",
    options=all_genes,
    default=[default_seed] if default_seed else []
)

# Baseline Metrics Calculation
base_nodes = base_graph.number_of_nodes()
base_edges = base_graph.number_of_edges()
base_cc = nx.number_connected_components(base_graph)
base_eff = nx.global_efficiency(base_graph)

st.subheader("1. Baseline Network Properties")
b1, b2, b3, b4 = st.columns(4)
b1.metric("Total Genes (Nodes)", base_nodes)
b2.metric("Total Edges", base_edges)
b3.metric("Connected Components", base_cc)
b4.metric("Global Efficiency", f"{base_eff:.4f}")

st.divider()

# Knockout Simulation Logic
if target_genes:
    st.subheader(f"2. Knockout Simulation Output: {', '.join(target_genes)}")
    
    # Perform node deletion
    perturbed_graph = base_graph.copy()
    perturbed_graph.remove_nodes_from(target_genes)
    
    post_nodes = perturbed_graph.number_of_nodes()
    post_edges = perturbed_graph.number_of_edges()
    edges_lost = base_edges - post_edges
    pct_edges_lost = (edges_lost / base_edges) * 100 if base_edges > 0 else 0
    
    post_cc = nx.number_connected_components(perturbed_graph)
    post_eff = nx.global_efficiency(perturbed_graph) if post_nodes > 0 else 0.0
    eff_loss = ((base_eff - post_eff) / base_eff) * 100 if base_eff > 0 else 0
    
    # Metrics Row
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Edges Removed", edges_lost, delta=f"-{pct_edges_lost:.1f}%")
    m2.metric("Remaining Nodes", post_nodes, delta=f"-{len(target_genes)}")
    m3.metric("Connected Components", post_cc, delta=f"+{post_cc - base_cc}")
    m4.metric("Global Efficiency Loss", f"{post_eff:.4f}", delta=f"-{eff_loss:.2f}%", delta_color="inverse")
    
    st.divider()
    
    # Degree Centrality Tables
    st.subheader("3. Hub Centrality Shift Analysis")
    c1, c2 = st.columns(2)
    
    with c1:
        st.write("#### Baseline Hub Ranking (Degree Centrality)")
        base_hubs = sorted(nx.degree_centrality(base_graph).items(), key=lambda x: x[1], reverse=True)[:8]
        st.dataframe(pd.DataFrame(base_hubs, columns=['Gene Symbol', 'Centrality Score']), use_container_width=True)
        
    with c2:
        st.write("#### Post-Knockout Remaining Hubs")
        post_hubs = sorted(nx.degree_centrality(perturbed_graph).items(), key=lambda x: x[1], reverse=True)[:8]
        st.dataframe(pd.DataFrame(post_hubs, columns=['Gene Symbol', 'Centrality Score']), use_container_width=True)
        
    st.divider()
    
    # Interactive Graph View
    st.subheader("4. Interactive Graph Visualizations")
    v1, v2 = st.columns(2)
    with v1:
        st.write("#### Baseline Topology")
        components.html(render_pyvis_network(base_graph, highlight_genes=target_genes), height=470)
    with v2:
        st.write("#### Post-Knockout Topology")
        components.html(render_pyvis_network(perturbed_graph), height=470)

else:
    st.info("Select one or more genes in the sidebar to run the structural removal test.")