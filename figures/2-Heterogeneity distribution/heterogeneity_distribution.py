import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

result_df = pd.read_table('/content/drive/MyDrive/GWAS/0_REVISED_VERSION/METAANALYSIS1_HETERO.TBL')
result_df

# Plot distribution of I-square
sns.distplot(result_df['HetISq'], hist=True, kde=False,
             bins=int(180/5), color = 'blue',
             hist_kws={'edgecolor':'black'},
             norm_hist=True)

plt.xlabel("I-square")
plt.ylabel("Frequency")
plt.show()

# Plot distribution of I-square with taking log to the y-axis
sns.distplot(result_df['HetISq'], hist=True, kde=False,
             bins=int(180/5), color = 'blue',
             hist_kws={'edgecolor':'black'},
             norm_hist=True)

plt.yscale('log')
plt.xlabel("I-square")
plt.ylabel("log10(Frequency)")
plt.show()

# Measure the percentage of significant p-value of Cochran's Q Test
total_snps = result_df.shape[0]
hetero_snps = result_df[result_df['HetPVal'] < 0.05].shape[0]
print(hetero_snps/total_snps*100)
