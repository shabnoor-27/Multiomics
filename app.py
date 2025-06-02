import os
import tempfile
import requests
import numpy as np
import pandas as pd
import streamlit as st
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px

from pyvis.network import Network
import networkx as nx
from gseapy import enrichr
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import umap
from matplotlib.colors import Normalize
import streamlit.components.v1 as components

# -----------------------------
# App Configuration
# -----------------------------
st.image("logo.png", width=200)
st.title("🧬 Multi-Omics Integration Vizzhy App")

with st.sidebar:
    st.markdown("---")
    st.markdown("**👨‍💻 Created by: SHABNOOR**")
    st.markdown("[LinkedIn](https://www.linkedin.com/in/priyadarshini24) | [GitHub](https://github.com/shabnoor-27)")


# -----------------------------
# File Upload Section
# -----------------------------
st.header("📁 Upload Omics Data")

genomics = st.file_uploader("Upload Genomics CSV", type="csv")
transcriptomics = st.file_uploader("Upload Transcriptomics CSV", type="csv")
proteomics = st.file_uploader("Upload Proteomics CSV", type="csv")


if genomics:
    gdf = pd.read_csv(genomics)

if transcriptomics:
    tdf = pd.read_csv(transcriptomics)

if proteomics:
    pdf = pd.read_csv(proteomics)

# -----------------------------
# Sidebar Filters
# -----------------------------
st.sidebar.header("⚙️ Settings")

log2FC_thresh = float(st.sidebar.text_input("Min log2FC Score (Genomics)", value="0.05"))
t_pVal_thresh = float(st.sidebar.text_input("Max p-value (Transcriptomics)", value="0.05"))
p_intensity_thresh = float(st.sidebar.text_input("Min Intensity (Proteomics)", value="100000"))


run_enrichment = st.sidebar.checkbox("Run Enrichment Analyses", value=True)
show_network = st.sidebar.checkbox("Show Network Visualization", value=True)
show_association_table = st.sidebar.checkbox("Show Association Table", value=True)
num_pathways_to_show = st.sidebar.slider("Number of Pathways to Display in Network", min_value=1, max_value=100, value=10)

# -----------------------------
# Preview Filtered Data: Top N Rows
# -----------------------------
preview_n = st.sidebar.slider("Preview Top N Filtered Rows", 5, 50, 10)

# Display Filtered Data Previews
st.subheader("🔍 Filtered Data Preview")

if genomics and transcriptomics and proteomics:
    try:
        gdf_filtered = gdf[gdf['log2FC'] >= log2FC_thresh]
        tdf_filtered = tdf[(tdf['pVal'] <= t_pVal_thresh)]
        
    # Filter Transcriptomics data based on logFC threshold
        if log2FC_thresh > 0.05:
    # For positive logFC threshold, filter for values >= logFC_thresh
              tdf_filtered = tdf[tdf['log2FC'] >= log2FC_thresh]
        elif log2FC_thresh < 0:
    # For negative logFC threshold, filter for values < logFC_thresh
              tdf_filtered = tdf[tdf['log2FC'] < log2FC_thresh]
        else:
    # If logFC_thresh is 0, no filtering needed
              tdf_filtered = tdf
        pdf_filtered = pdf[pdf['Intensity'] >= p_intensity_thresh]
        
        # Display filtered data for all three omics
        st.markdown("**Genomics**")
        st.dataframe(gdf_filtered.head(preview_n))  # Preview the top N filtered rows for genomics
        st.markdown("**Transcriptomics**")
        st.dataframe(tdf_filtered.head(preview_n))  # Preview the top N filtered rows for transcriptomics
        st.markdown("**Proteomics**")
        st.dataframe(pdf_filtered.head(preview_n))  # Preview the top N filtered rows for proteomics

    except Exception as e:
        st.error(f"Integration error: {e}")

# -----------------------------
# Filtering and Integration
# -----------------------------
st.header("🎛️ Filter & Integrate")

if genomics and transcriptomics and proteomics:
    try:
        gdf_filtered = gdf[gdf['log2FC'] >= log2FC_thresh]
        tdf_filtered = tdf[(tdf['pVal'] <= t_pVal_thresh)]
        pdf_filtered = pdf[pdf['Intensity'] >= p_intensity_thresh]

        union_genes = set(gdf_filtered['Gene']) | set(tdf_filtered['Gene'])

        def extract_uniprot_ids(protein_series):
            ids = set()
            for entry in protein_series.dropna():
                for pid in str(entry).split(";"):
                    if pid.strip():
                        ids.add(pid.strip())
            return ids

        def map_uniprot_to_gene(uniprot_ids):
            mapping = {}
            ids = list(uniprot_ids)
            for i in range(0, len(ids), 100):
                chunk = ids[i:i+100]
                query = " OR ".join([f"accession:{id_}" for id_ in chunk])
                url = f"https://rest.uniprot.org/uniprotkb/search?query={query}&fields=accession,gene_names&format=tsv"
                try:
                    r = requests.get(url)
                    if r.status_code == 200:
                        lines = r.text.strip().split('\n')[1:]
                        for line in lines:
                            acc, genes = line.split('\t')
                            mapping[acc] = genes.split()[0] if genes else acc
                except Exception as e:
                    st.warning(f"UniProt API error: {e}")
            return mapping

        unique_uniprot_ids = extract_uniprot_ids(pdf_filtered['Protein IDs'])
        uniprot_gene_map = map_uniprot_to_gene(unique_uniprot_ids)

        expanded_rows = []
        for _, row in pdf_filtered.iterrows():
            for pid in str(row['Protein IDs']).split(';'):
                pid = pid.strip()
                gene = uniprot_gene_map.get(pid)
                if gene:
                    expanded_rows.append({'Protein IDs': pid, 'GeneName': gene})

        expanded_protein_df = pd.DataFrame(expanded_rows)
        protein_gene_map = dict(zip(expanded_protein_df['Protein IDs'], expanded_protein_df['GeneName']))

        all_entities = union_genes | set(protein_gene_map.values())

        results = {}
        raw_assoc_data = []

        if run_enrichment:
            st.header("📊 Enrichment Analyses")
            libraries = {
                "Mouse Disease Ontology": "DO_Mouse_2021",
                "Disease Associations": "MGI_Mammalian_Phenotype_Level_4_2021",
                "KEGG Mouse Metabolites": "KEGG_2019_Mouse",
                "KEGG Mouse Metabolites": "KEGG_2019_Mouse",
                "WikiPathways (Mouse)": "WikiPathways_2019_Mouse"
            }


            for name, lib in libraries.items():
                try:
                    gene_list_clean = [str(g).strip() for g in union_genes if pd.notna(g)]
                    enr = enrichr(gene_list=gene_list_clean, gene_sets=lib, outdir=None)

                    if enr.results.empty:
                        continue

                    df = enr.results.copy()
                    df['-log10(pval)'] = -np.log10(df['P-value'])
                    df = df.rename(columns={"Term": "Pathway", "Genes": "Genes_Involved"})
                    results[name] = df

                    fig = px.bar(df.head(10), x="Pathway", y="-log10(pval)", title=f"Top 10 {name}")
                    st.plotly_chart(fig)
                except Exception as e:
                    st.error(f"Error in {name}: {e}")

        if show_network and results:
            st.subheader("🧠 Interactive Omics Network")
            net = Network(height='800px', width='100%', directed=False)
            net.force_atlas_2based()

            legend_items = {
                "Gene": 'gray-+', "Protein IDs": 'gold',
                "Pathway": 'skyblue', "Metabolite": 'lightgreen', "Disease": 'lightcoral'
            }

            for i, (label, color) in enumerate(legend_items.items()):
                net.add_node(f"legend_{label}", label=label, shape='box', color=color, size=20, x=-1000, y=-i*50, physics=False, fixed=True)

            color_map = {
                "Mouse Disease Ontology": "skyblue",
                "Disease Associations": "lightcoral",
                "KEGG Mouse Metabolites": "pink",
                "KEGG Mouse Metabolites": "orange",
                "WikiPathways (Mouse)": "lightgreen"
            }

            for name, df in results.items():
                color = color_map.get(name, "gray")
                for _, row in df.head(num_pathways_to_show).iterrows():
                    term = row['Pathway']
                    net.add_node(term, label=term, color=color)
                    for gene in row['Genes_Involved'].split(';'):
                        gene = gene.strip()
                        if not gene:
                            continue
                        net.add_node(gene, label=gene, color='gray')
                        net.add_edge(gene, term)
                        matched_proteins = [prot for prot, gname in protein_gene_map.items() if gname == gene]
                        for prot in matched_proteins:
                            net.add_node(prot, label=prot, color='gold')
                            net.add_edge(gene, prot)
                        raw_assoc_data.append({
                            'Gene': gene,
                            'Protein IDs': ';'.join(matched_proteins),
                            'Pathway': term if name == 'Reactome Pathways' else '',
                            'Metabolite': term if name == 'HMDB Metabolites' else '',
                            'Disease': term if name == 'Disease Associations' else ''
                        })

            with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as tmp_file:
                net.save_graph(tmp_file.name)
                st.components.v1.html(open(tmp_file.name, 'r', encoding='utf-8').read(), height=800)

       if show_association_table and raw_assoc_data:
        st.subheader("📄 Gene-Protein-Term Association Summary")

        # Convert raw data to DataFrame
        df = pd.DataFrame(raw_assoc_data)

        # Ensure required columns exist
        required_columns = ['Gene', 'Protein IDs', 'Pathway', 'Disease', 'Metabolite']
        if not all(col in df.columns for col in required_columns):
            st.warning("Some required columns are missing in the data.")
            st.write("Expected columns:", required_columns)
            st.write("Found columns:", list(df.columns))
        else:
            # Group and aggregate data
            assoc_df = df.groupby('Gene').agg({
                'Protein IDs': lambda x: ';'.join(set(filter(None, map(str, x)))),
                'Pathway': lambda x: ';'.join(set(filter(None, map(str, x)))),
                'Disease': lambda x: ';'.join(set(filter(None, map(str, x)))),
                'Metabolite': lambda x: ';'.join(set(filter(None, map(str, x))))
            }).reset_index()

            # Count non-null fields to sort
            assoc_df['non_nulls'] = assoc_df.replace('', pd.NA).notnull().sum(axis=1)
            assoc_df = assoc_df.sort_values(by='non_nulls', ascending=False).drop(columns='non_nulls')

            # Display table
            st.dataframe(assoc_df)

except Exception as e:
    st.error(f"Integration error: {e}")
# -----------------------------
# UMAP + KMeans Clustering
# -----------------------------
# -----------------------------
# UMAP + KMeans Clustering
# -----------------------------
st.header("📉 UMAP Clustering")

if genomics and transcriptomics and proteomics:
    try:
        # Make sure filtered data exists
        if 'gdf_filtered' in locals() and 'tdf_filtered' in locals() and 'pdf_filtered' in locals():
            # Merge datasets on Gene column
            merged_df = pd.merge(gdf_filtered[['Gene', 'log2FC']], tdf_filtered[['Gene', 'log2FC']], on='Gene', suffixes=('_genomics', '_transcriptomics'))
            merged_df = pd.merge(merged_df, pdf_filtered[['Gene', 'Intensity']], on='Gene')

            # Prepare data for UMAP
            features = merged_df[['log2FC_genomics', 'log2FC_transcriptomics', 'Intensity']]
            features_scaled = StandardScaler().fit_transform(features)

            reducer = umap.UMAP(random_state=42)
            embedding = reducer.fit_transform(features_scaled)

            kmeans = KMeans(n_clusters=3, random_state=42).fit(embedding)
            merged_df['Cluster'] = kmeans.labels_

            fig = px.scatter(x=embedding[:, 0], y=embedding[:, 1], color=merged_df['Cluster'].astype(str),
                             labels={'x': 'UMAP 1', 'y': 'UMAP 2'}, title="UMAP Clustering of Integrated Omics Data")
            st.plotly_chart(fig)
        else:
            st.warning("Filtered data not found. Please make sure the data integration step ran without error.")

    except Exception as e:
        st.error(f"UMAP clustering error: {e}")
