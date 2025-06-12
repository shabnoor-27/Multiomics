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
st.markdown("**👨‍💻 Created by: SHABNOOR**")
st.markdown("[LinkedIn](https://www.linkedin.com/in/priyadarshini24) | [GitHub](https://github.com/shabnoor-27)")
 
 
 
# -----------------------------
 
# File Upload Section
 
# -----------------------------
 
st.header("📁 Upload Omics Data")
 
 
 
genomics = st.file_uploader("Upload Genomics CSV", type="csv")
 
transcriptomics = st.file_uploader("Upload Transcriptomics CSV", type="csv")
 
proteomics = st.file_uploader("Upload Proteomics CSV", type="csv")
 
 
 
gdf = tdf = pdf = None
 
 
 
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
 
 
 
st.subheader("🔍 Filtered Data Preview")
 
 
 
if gdf is not None and tdf is not None and pdf is not None:
 
try:
 
# Filter Genomics
 
gdf_filtered = gdf[gdf['log2FC'] >= log2FC_thresh]
 
 
 
# Filter Transcriptomics
 
if log2FC_thresh > 0.05:
 
tdf_filtered = tdf[(tdf['pVal'] <= t_pVal_thresh) & (tdf['log2FC'] >= log2FC_thresh)]
 
elif log2FC_thresh < 0:
 
tdf_filtered = tdf[(tdf['pVal'] <= t_pVal_thresh) & (tdf['log2FC'] < log2FC_thresh)]
 
else:
 
tdf_filtered = tdf[tdf['pVal'] <= t_pVal_thresh]
 
 
 
# Handle Proteomics: check for protein ID column
 
prot_id_col = None
 
for col_candidate in ['Protein IDs', 'Protein', 'ProteinID', 'Protein_Id']:
 
if col_candidate in pdf.columns:
 
prot_id_col = col_candidate
 
break
 
 
 
if prot_id_col is None:
 
st.error("Could not find a protein ID column in proteomics data. Please ensure the column is named 'Protein IDs' or similar.")
 
st.stop()
 
 
 
# Filter Proteomics by Intensity and require Gene column
 
if 'Gene' not in pdf.columns:
 
st.error("Proteomics data must contain a 'Gene' column to match proteins with genes.")
 
st.stop()
 
 
 
pdf_filtered = pdf[pdf['Intensity'] >= p_intensity_thresh]
 
 
 
# Show filtered data previews
 
st.markdown("**Genomics**")
 
st.dataframe(gdf_filtered.head(preview_n))
 
 
 
st.markdown("**Transcriptomics**")
 
st.dataframe(tdf_filtered.head(preview_n))
 
 
 
st.markdown("**Proteomics**")
 
st.dataframe(pdf_filtered.head(preview_n))
 
 
 
except Exception as e:
 
st.error(f"Integration error: {e}")
 
 
 
# -----------------------------
 
# Filtering and Integration
 
# -----------------------------
 
st.header("🎛️ Filter & Integrate")
 
 
 
if gdf is not None and tdf is not None and pdf is not None:
 
try:
 
# Repeat filtering for integration
 
gdf_filtered = gdf[gdf['log2FC'] >= log2FC_thresh]
 
 
 
if log2FC_thresh > 0.05:
 
tdf_filtered = tdf[(tdf['pVal'] <= t_pVal_thresh) & (tdf['log2FC'] >= log2FC_thresh)]
 
elif log2FC_thresh < 0:
 
tdf_filtered = tdf[(tdf['pVal'] <= t_pVal_thresh) & (tdf['log2FC'] < log2FC_thresh)]
 
else:
 
tdf_filtered = tdf[tdf['pVal'] <= t_pVal_thresh]
 
 
 
pdf_filtered = pdf[pdf['Intensity'] >= p_intensity_thresh]
 
 
 
# Get union of genes
 
union_genes = set(gdf_filtered['Gene']) | set(tdf_filtered['Gene'])
 
 
 
# Extract protein IDs from proteomics
 
prot_id_col = None
 
for col_candidate in ['Protein IDs', 'Protein', 'ProteinID', 'Protein_Id']:
 
if col_candidate in pdf.columns:
 
prot_id_col = col_candidate
 
break
 
 
 
if prot_id_col is None:
 
st.error("Could not find a protein ID column in proteomics data. Please ensure the column is named 'Protein IDs' or similar.")
 
st.stop()
 
 
 
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
 
 
 
unique_uniprot_ids = extract_uniprot_ids(pdf_filtered[prot_id_col])
 
uniprot_gene_map = map_uniprot_to_gene(unique_uniprot_ids)
 
 
 
expanded_rows = []
 
for _, row in pdf_filtered.iterrows():
 
for pid in str(row[prot_id_col]).split(';'):
 
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
 
"Gene": 'gray', "Protein IDs": 'gold',
 
"Pathway": 'skyblue', "Metabolite": 'lightgreen', "Disease": 'lightcoral'
 
}
 
 
 
for i, (label, color) in enumerate(legend_items.items()):
 
net.add_node(f"legend_{label}", label=label, shape='box', color=color, size=20, x=-1000, y=-i*50, physics=False, fixed=True)
 
 
 
color_map = {
 
"Mouse Disease Ontology": "skyblue",
 
"Disease Associations": "lightcoral",
 
"KEGG Mouse Metabolites": "lightgreen",
 
"WikiPathways (Mouse)": "orange"
 
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
 
'Pathway': term if name == 'Mouse Disease Ontology' else '',
 
'Disease': term if name == 'Disease Associations' else '',
 
'Metabolite': term if name == 'KEGG Mouse Metabolites' else ''
 
})
 
 
 
with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as tmp_file:
 
net.save_graph(tmp_file.name)
 
html = open(tmp_file.name, 'r', encoding='utf-8').read()
 
st.components.v1.html(html, height=800)
 
 
 
if show_association_table and raw_assoc_data:
 
st.subheader("📄 Gene-Protein-Term Association Summary")
 
 
 
df_assoc = pd.DataFrame(raw_assoc_data)
 
 
 
required_columns = ['Gene', 'Protein IDs', 'Pathway', 'Disease', 'Metabolite']
 
if not all(col in df_assoc.columns for col in required_columns):
 
st.warning("Some required columns are missing in the association data.")
 
else:
 
assoc_df = df_assoc.groupby('Gene').agg({
 
'Protein IDs': lambda x: ';'.join(set(filter(None, map(str, x)))),
 
'Pathway': lambda x: ';'.join(set(filter(None, map(str, x)))),
 
'Disease': lambda x: ';'.join(set(filter(None, map(str, x)))),
 
'Metabolite': lambda x: ';'.join(set(filter(None, map(str, x))))
 
}).reset_index()
 
 
 
assoc_df['non_nulls'] = assoc_df.replace('', pd.NA).notnull().sum(axis=1)
 
assoc_df = assoc_df.sort_values(by='non_nulls', ascending=False).drop(columns='non_nulls')
 
 
 
st.dataframe(assoc_df)
 
 
 
except Exception as e:
 
st.error(f"Integration error: {e}")
 
 
 
# -----------------------------
 
# UMAP + KMeans Clustering
 
# -----------------------------
 
st.header("📉 UMAP Clustering")
 
 
 
if gdf is not None and tdf is not None and pdf is not None:
 
try:
 
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
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
2. code:  
 
 
 
 
 
 
 
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
 
 
 
gdf = tdf = pdf = None
 
 
 
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
 
 
 
st.subheader("🔍 Filtered Data Preview")
 
 
 
if gdf is not None and tdf is not None and pdf is not None:
 
try:
 
gdf_filtered = gdf[gdf['log2FC'] >= log2FC_thresh]
 
 
 
if log2FC_thresh > 0.05:
 
tdf_filtered = tdf[(tdf['pVal'] <= t_pVal_thresh) & (tdf['log2FC'] >= log2FC_thresh)]
 
elif log2FC_thresh < 0:
 
tdf_filtered = tdf[(tdf['pVal'] <= t_pVal_thresh) & (tdf['log2FC'] < log2FC_thresh)]
 
else:
 
tdf_filtered = tdf[tdf['pVal'] <= t_pVal_thresh]
 
 
 
prot_id_col = None
 
for col_candidate in ['Protein IDs', 'Protein', 'ProteinID', 'Protein_Id']:
 
if col_candidate in pdf.columns:
 
prot_id_col = col_candidate
 
break
 
 
 
if prot_id_col is None:
 
st.error("Could not find a protein ID column in proteomics data. Please ensure the column is named 'Protein IDs' or similar.")
 
st.stop()
 
 
 
if 'Gene' not in pdf.columns:
 
st.error("Proteomics data must contain a 'Gene' column to match proteins with genes.")
 
st.stop()
 
 
 
pdf_filtered = pdf[pdf['Intensity'] >= p_intensity_thresh]
 
 
 
st.markdown("**Genomics**")
 
st.dataframe(gdf_filtered.head(preview_n))
 
 
 
st.markdown("**Transcriptomics**")
 
st.dataframe(tdf_filtered.head(preview_n))
 
 
 
st.markdown("**Proteomics**")
 
st.dataframe(pdf_filtered.head(preview_n))
 
 
 
except Exception as e:
 
st.error(f"Integration error: {e}")
 
 
 
# -----------------------------
 
# Filtering and Integration
 
# -----------------------------
 
st.header("🎛️ Filter & Integrate")
 
 
 
if gdf is not None and tdf is not None and pdf is not None:
 
try:
 
gdf_filtered = gdf[gdf['log2FC'] >= log2FC_thresh]
 
 
 
if log2FC_thresh > 0.05:
 
tdf_filtered = tdf[(tdf['pVal'] <= t_pVal_thresh) & (tdf['log2FC'] >= log2FC_thresh)]
 
elif log2FC_thresh < 0:
 
tdf_filtered = tdf[(tdf['pVal'] <= t_pVal_thresh) & (tdf['log2FC'] < log2FC_thresh)]
 
else:
 
tdf_filtered = tdf[tdf['pVal'] <= t_pVal_thresh]
 
 
 
pdf_filtered = pdf[pdf['Intensity'] >= p_intensity_thresh]
 
 
 
union_genes = set(gdf_filtered['Gene']) | set(tdf_filtered['Gene'])
 
 
 
prot_id_col = None
 
for col_candidate in ['Protein IDs', 'Protein', 'ProteinID', 'Protein_Id']:
 
if col_candidate in pdf.columns:
 
prot_id_col = col_candidate
 
break
 
 
 
if prot_id_col is None:
 
st.error("Could not find a protein ID column in proteomics data. Please ensure the column is named 'Protein IDs' or similar.")
 
st.stop()
 
 
 
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
 
 
 
unique_uniprot_ids = extract_uniprot_ids(pdf_filtered[prot_id_col])
 
uniprot_gene_map = map_uniprot_to_gene(unique_uniprot_ids)
 
 
 
expanded_rows = []
 
for _, row in pdf_filtered.iterrows():
 
for pid in str(row[prot_id_col]).split(';'):
 
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
 
"Gene": 'gray', "Protein IDs": 'gold',
 
"Pathway": 'skyblue', "Metabolite": 'lightgreen', "Disease": 'lightcoral'
 
}
 
 
 
for i, (label, color) in enumerate(legend_items.items()):
 
net.add_node(f"legend_{label}", label=label, shape='box', color=color, size=20, x=-1000, y=-i*50, physics=False, fixed=True)
 
 
 
color_map = {
 
"Mouse Disease Ontology": "skyblue",
 
"Disease Associations": "lightcoral",
 
"KEGG Mouse Metabolites": "lightgreen",
 
"WikiPathways (Mouse)": "orange"
 
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
 
'Pathway': term if name == 'Mouse Disease Ontology' else '',
 
'Disease': term if name == 'Disease Associations' else '',
 
'Metabolite': term if name == 'KEGG Mouse Metabolites' else ''
 
})
 
 
 
with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as tmp_file:
 
net.save_graph(tmp_file.name)
 
html = open(tmp_file.name, 'r', encoding='utf-8').read()
 
st.components.v1.html(html, height=800)
 
 
 
if show_association_table and raw_assoc_data:
 
st.subheader("📄 Gene-Protein-Term Association Summary")
 
 
 
df_assoc = pd.DataFrame(raw_assoc_data)
 
 
 
for col in ['Gene', 'Protein IDs', 'Pathway', 'Disease', 'Metabolite']:
 
if col not in df_assoc.columns:
 
df_assoc[col] = ""
 
 
 
assoc_df = df_assoc.groupby('Gene').agg({
 
'Protein IDs': lambda x: ';'.join(set(pid for pid in x if pid)),
 
'Pathway': lambda x: ';'.join(set(p for p in x if p)),
 
'Disease': lambda x: ';'.join(set(d for d in x if d)),
 
'Metabolite': lambda x: ';'.join(set(m for m in x if m))
 
}).reset_index()
 
 
 
assoc_df['non_nulls'] = assoc_df.replace('', pd.NA).notnull().sum(axis=1)
 
assoc_df = assoc_df.sort_values(by='non_nulls', ascending=False).drop(columns='non_nulls')
 
 
 
st.dataframe(assoc_df)
 
 
 
except Exception as e:
 
st.error(f"Integration error: {e}")
 
 
 
# -----------------------------
 
# UMAP + KMeans Clustering
 
# -----------------------------
 
st.header("📉 UMAP Clustering")
 
 
 
if gdf is not None and tdf is not None and pdf is not None:
 
try:
 
if 'gdf_filtered' in locals() and 'tdf_filtered' in locals() and 'pdf_filtered' in locals():
 
merged_df = pd.merge(gdf_filtered[['Gene', 'log2FC']], tdf_filtered[['Gene', 'log2FC']], on='Gene', suffixes=('_genomics', '_transcriptomics'))
 
merged_df = pd.merge(merged_df, pdf_filtered[['Gene', 'Intensity']], on='Gene')
 
 
 
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
