import pandas as pd
from Ensembl_converter import EnsemblConverter
from scipy import stats

# Load TWAS results and clean columns
twas_df = pd.read_table("/content/drive/MyDrive/GWAS/TWAS/twas_result.tsv.GW").drop(columns=['PANEL', 'FILE']).dropna()
twas_df['TWAS.P'] = twas_df['TWAS.P'].astype(float)
twas_df['Gene_ID'] = twas_df['ID'].str.split('.').str[0]

# Convert Ensembl IDs to gene symbols
converter = EnsemblConverter(use_progress_bar=True)
result = converter.convert_ids(twas_df['Gene_ID'])

# Merge conversion results
preprocessed_twas_df = result.merge(twas_df, left_on='ENSG', right_on='Gene_ID').drop(columns=['ID', 'Gene_ID'])

# Sort by CHR and P0
preprocessed_twas_df = preprocessed_twas_df.sort_values(['CHR', 'P0'])

# FDR correction
preprocessed_twas_df['TWAS.FDR'] = stats.false_discovery_control(preprocessed_twas_df['TWAS.P'])

# Save preprocessed data and significant results
preprocessed_twas_df.to_csv("/content/drive/MyDrive/GWAS/TWAS/preprocessed_twas.csv", index=False)
preprocessed_twas_df.query('TWAS.FDR < 0.05')[['ENSG', 'Symbol']].to_csv("/content/drive/MyDrive/GWAS/TWAS/sig_fdr_twas.txt", index=False)