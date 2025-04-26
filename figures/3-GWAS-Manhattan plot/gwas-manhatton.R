library(dplyr)
library(qqman)
library(RColorBrewer)

# Load and save GWAS data
gwasResults <- read.table("METAANALYSIS1.TBL", header=TRUE) %>%
  rename(BP = POS)
write.csv(gwasResults, "METAANALYSIS1.csv", row.names=FALSE)

# Chromosome summary
table(gwasResults$CHR)

# Manhattan plot
tiff("manhattan_plot.tiff", units="in", width=30, height=15, res=700)
manhattan(gwasResults, col=brewer.pal(8, "Set2"), cex=2, cex.axis=1, ylim=c(0, 35),
          genomewideline = -log10(5e-8), suggestiveline = -log10(1e-5))
dev.off()
