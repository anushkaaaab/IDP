# Gene Regulatory Network (GRN) Modeling & Dynamic Perturbation Analysis in Breast Cancer

An end-to-end bioinformatics and network biology framework designed to construct breast cancer gene regulatory networks from transcriptomic profiles, calculate node centrality metrics, and simulate *in silico* knockouts to observe network topological cascades.

---

##  Project Overview

Targeting key regulator genes in cancer therapy often causes unpredictable system-wide network shifts. This project addresses that challenge by building a functional bioinformatics pipeline that:
1. **Ingests & Normalizes** transcriptomic expression data from public cohorts (NCBI GEO).
2. **Infers Gene Regulatory Networks (GRNs)** using pairwise correlation and mutual information thresholding.
3. **Identifies Key Hub Genes** using graph centrality metrics (degree, betweenness, eigenvector).
4. **Simulates *In Silico* Knockouts** to analyze structural stability and signal flow disruption.
5. **Provides an Interactive Dashboard** for non-technical researchers to explore network dynamics.

---

##  Project Status (Review 2 Milestone — 20% Complete)

| Module / Milestone | Description | Status |
| :--- | :--- | :--- |
| **Stage 1: Data Ingestion & Mapping** | Automated pipeline via `GEOparse`, probe-to-HUGO symbol mapping, log2 RMA normalization |  **Completed** |
| **Stage 1: Expression Processing** | Missing value handling, duplicate probe aggregation, variance filtering |  **Completed** |
| **Stage 2: Baseline GRN Inference** | Pearson co-expression matrix calculation, sparsification ($\vert{}r\vert{} \ge 0.65$), NetworkX graph generation |  **Completed** |
| **Stage 2: Centrality Analysis** | Initial degree centrality scoring & top hub gene identification |  **Completed** |
| **Stage 3: Perturbation Engine** | Structural removal logic (*in silico* gene knockout simulation) | ⏳ In Progress |
| **Stage 4: Interactive UI** | Streamlit/PyVis interactive graph explorer & dashboard | ⏳ In Progress |

---

##  Dataset & Platform

* **Dataset:** [GSE3494](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3494) (Uppsala Breast Cancer Cohort)
* **Platform:** Affymetrix Human Genome U133 Array Set (`GPL96` / `GPL97`)
* **Cohort Size:** $N = 251$ primary breast cancer tissue samples
* **Target Driver Genes Analyzed:** *TP53*, *ESR1*, *ERBB2*, *BRCA1*, *BRCA2*, *PGR*

---

## 🛠️ System Architecture
