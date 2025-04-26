import pandas as pd

# Load files
metal_df = pd.read_table('/content/drive/MyDrive/GWAS/METAL/METAANALYSIS1.TBL')
gwas_catalog_df = pd.read_table('/content/drive/MyDrive/GWAS/METAL/EFO_0004214_associations_export.tsv')
ind_sig_snp_df = pd.read_table('/content/drive/MyDrive/GWAS/FUMA/SNP2GENE/Result Files/IndSigSNPs.txt')
mapped_gene_df = pd.read_table('/content/drive/MyDrive/GWAS/FUMA/SNP2GENE/Result Files/genes.txt')
genecards_aaa_df = pd.read_csv('/content/drive/MyDrive/GWAS/METAL/GeneCards-SearchResults.csv')

# Define thresholds
genome_wide_sig, suggestive_sig = 5e-8, 1e-5

# Filter significant SNPs
sig1_metal_df = metal_df[metal_df['P-value'] < genome_wide_sig]
sig2_metal_df = metal_df[metal_df['P-value'] < suggestive_sig]

# Identify novel SNPs
snp_catalog = set(gwas_catalog_df['riskAllele'].str[:-2])
novel_snp = set(ind_sig_snp_df['rsID']) - snp_catalog
novel_snp_df = ind_sig_snp_df[ind_sig_snp_df['rsID'].isin(novel_snp)]
novel_snp_df.to_csv('/content/drive/MyDrive/GWAS/METAL/novel_snp.csv', index=False)

# Identify novel genes
catalog_genes = set(j for i in gwas_catalog_df['mappedGenes'].dropna() for j in i.split(','))
fuma_genes = set(mapped_gene_df['symbol'])
genecards_genes = set(genecards_aaa_df['Gene Symbol'])

novel_genes = fuma_genes - catalog_genes - genecards_genes
novel_gene_df = mapped_gene_df[mapped_gene_df['symbol'].isin(novel_genes)]
novel_gene_df.to_csv('/content/drive/MyDrive/GWAS/METAL/novel_gene.csv', index=False)