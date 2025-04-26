# Import packages
import pandas as pd
import numpy as np
import networkx as nx

# Load DEGs from TWAS and FUMA
sig_twas = pd.read_csv("/content/drive/MyDrive/GWAS/TWAS/sig_fdr_twas.txt")
eqtl = pd.read_table('/content/drive/MyDrive/GWAS/Network/eqtl.txt')

# Combine DEGs from TWAS and eQTLs
degs = pd.Series(pd.concat([sig_twas['Symbol'], eqtl['symbol']]).unique())
degs.to_csv('/content/drive/MyDrive/GWAS/Network/DEGs.txt', index=False, header=False)

# Download and read STRING interactome and protein info
!wget -q https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz
!gunzip -f /content/9606.protein.links.v12.0.txt.gz
ppi = pd.read_table("/content/9606.protein.links.v12.0.txt", sep=' ')

!wget -q https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz
!gunzip -f /content/9606.protein.info.v12.0.txt.gz
prot_info = pd.read_table("/content/9606.protein.info.v12.0.txt")[['#string_protein_id', 'preferred_name']]

# Map STRING IDs to gene symbols
ppi = ppi.merge(prot_info, left_on='protein2', right_on='#string_protein_id') \
         .rename(columns={'preferred_name': 'protein2'}).drop(columns=['#string_protein_id'])
ppi = ppi.merge(prot_info, left_on='protein1', right_on='#string_protein_id') \
         .rename(columns={'preferred_name': 'protein1'}).drop(columns=['#string_protein_id'])

# Filter PPI network to DEGs
deg_list = pd.read_table('/content/drive/MyDrive/GWAS/Network/DEGs.txt', header=None)[0]
ppi = ppi[ppi['protein1'].isin(deg_list) & ppi['protein2'].isin(deg_list)]

# Add weight columns
ppi['relationship_weight'] = ppi['combined_score'] / 1000
ppi['distance_weight'] = 1 - ppi['relationship_weight']

# Create graphs
G_rel = nx.from_pandas_edgelist(ppi, 'protein1', 'protein2', edge_attr='relationship_weight')
G_dist = nx.from_pandas_edgelist(ppi, 'protein1', 'protein2', edge_attr='distance_weight')

# Compute centralities
degree = {k: v/max(dict(G_rel.degree(weight='relationship_weight')).values()) for k, v in G_rel.degree(weight='relationship_weight')}
betweenness = nx.betweenness_centrality(G_dist, normalized=True, weight='distance_weight')

# Merge centralities
central_df = pd.DataFrame({'protein': degree.keys(), 'degree': degree.values()}).merge(
    pd.DataFrame({'protein': betweenness.keys(), 'betweenness': betweenness.values()}), on='protein')

# Save centrality results
central_df.to_csv('/content/drive/MyDrive/GWAS/Network/centrality.csv', index=False)

# 90th percentile proteins
p90_deg = central_df[central_df['degree'] >= np.percentile(central_df['degree'], 90)]
p90_btw = central_df[central_df['betweenness'] >= np.percentile(central_df['betweenness'], 90)]
