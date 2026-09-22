import streamlit as st
import pandas as pd
import networkx as nx

st.set_page_config(page_title="GRN Cancer Gene Perturbation Tool", layout="wide")

st.title("Gene Regulatory Network (GRN) Perturbation Simulator")
st.markdown("### Breast Cancer (GSE3494) Co-Expression Network Analysis")

@st.cache_data
def load_graph(csv_path):
    df = pd.read_csv(csv_path)
    G = nx.Graph()
    for _, row in df.iterrows():
        G.add_edge(row['gene_1'], row['gene_2'], weight=row['weight'])
    return G

try:
    base_graph = load_graph("network_edges.csv")
except Exception as e:
    st.error("Error loading 'network_edges.csv'. Ensure Student A's code has run.")
    st.stop()

st.sidebar.header("Perturbation Parameters")
all_genes = sorted(list(base_graph.nodes()))
target_genes = st.sidebar.multiselect(
    "Select Target Gene(s) to Knock Out:",
    options=all_genes,
    default=[all_genes[0]] if all_genes else []
)

# Baseline Summary Metrics
col1, col2, col3 = st.columns(3)
with col1:
    st.metric("Total Genes (Nodes)", base_graph.number_of_nodes())
with col2:
    st.metric("Total Interactions (Edges)", base_graph.number_of_edges())
with col3:
    st.metric("Connected Components", nx.number_connected_components(base_graph))

st.divider()

# Perturbation Engine Logic Execution
if target_genes:
    st.subheader(f"Knockout Simulation: {', '.join(target_genes)}")
    
    # Execute knockout on graph copy
    perturbed_graph = base_graph.copy()
    perturbed_graph.remove_nodes_from(target_genes)
    
    # Calculate Disruption Metrics
    orig_edges = base_graph.number_of_edges()
    new_edges = perturbed_graph.number_of_edges()
    edges_lost = orig_edges - new_edges
    pct_lost = (edges_lost / orig_edges) * 100 if orig_edges > 0 else 0
    
    orig_cc = nx.number_connected_components(base_graph)
    new_cc = nx.number_connected_components(perturbed_graph)

    # Display Metrics Side-by-Side
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Edges Removed", edges_lost, delta=f"-{pct_lost:.1f}%")
    m2.metric("Remaining Nodes", perturbed_graph.number_of_nodes(), delta=f"-{len(target_genes)}")
    m3.metric("New Components", new_cc, delta=f"+{new_cc - orig_cc}")
    
    # Identify Top Hubs in Perturbed Graph
    deg_centrality = nx.degree_centrality(perturbed_graph)
    top_remaining_hubs = sorted(deg_centrality.items(), key=lambda x: x[1], reverse=True)[:5]
    
    with m4:
        st.markdown("**Top Remaining Hubs:**")
        for gene, score in top_remaining_hubs:
            st.text(f"• {gene} ({score:.2f})")

    st.divider()

    # Network Comparison View
    st.subheader("Network Comparison")
    c1, c2 = st.columns(2)
    
    with c1:
        st.write("#### Baseline Network")
        st.info(f"Nodes: {base_graph.number_of_nodes()} | Edges: {base_graph.number_of_edges()}")
        # Display top 10 highest-degree nodes table
        base_hubs = sorted(nx.degree_centrality(base_graph).items(), key=lambda x: x[1], reverse=True)[:10]
        st.dataframe(pd.DataFrame(base_hubs, columns=['Gene', 'Degree Centrality']), use_container_width=True)
        
    with c2:
        st.write("#### Post-Knockout Network")
        st.warning(f"Nodes: {perturbed_graph.number_of_nodes()} | Edges: {perturbed_graph.number_of_edges()}")
        post_hubs = sorted(deg_centrality.items(), key=lambda x: x[1], reverse=True)[:10]
        st.dataframe(pd.DataFrame(post_hubs, columns=['Gene', 'Degree Centrality']), use_container_width=True)

else:
    st.info("Select one or more genes from the sidebar to execute the perturbation simulation.")