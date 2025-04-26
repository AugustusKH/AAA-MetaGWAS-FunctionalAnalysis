library(dplyr)
library(qqman)
library(RColorBrewer)

# Load and save GWAS data
gwasResults <- read.table("METAANALYSIS1.TBL", header=TRUE) %>%
  rename(BP = POS)
write.csv(gwasResults, "METAANALYSIS1.csv", row.names=FALSE)

# Chromosome summary
table(gwasResults$CHR)

# QQ plot
tiff("qq_plot.tiff", units="in", width=6, height=5, res=500)
qq(gwasResults$P)
dev.off()
