library(dplyr)
library(qqman)
library(RColorBrewer)

# Load and save GWAS data
twasResults <- read.table("twas_result.tsv.GW", header=TRUE) %>%
  rename(BP = POS)
write.csv(twasResults, "twas.csv", row.names=FALSE)

# Manhattan plot
tiff("manhattan_plot.tiff", units="in", width=30, height=15, res=700)
manhattan(twasResults, col=brewer.pal(8, "Set2"), cex=2, cex.axis=1, ylim=c(0, 35),
          genomewideline = -log10(5e-8), suggestiveline = -log10(1e-5))
dev.off()
