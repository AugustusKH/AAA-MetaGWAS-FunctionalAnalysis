library(dplyr)
library(qqman)
library(RColorBrewer)

# Load and save GWAS data
twasResults <- read.table("twas_result.tsv.GW", header=TRUE) %>%
  rename(BP = POS)
write.csv(twasResults, "twas.csv", row.names=FALSE)

# QQ plot
tiff("qq_plot.tiff", units="in", width=6, height=5, res=500)
qq(twasResults$P)
dev.off()
