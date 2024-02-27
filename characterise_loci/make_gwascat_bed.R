# Author: Samvida S. Venkatesh
# Date: 01/11/2021

library(tidyverse)

filepath <- "/well/lindgren/samvida/Resources/GWASCatalog"

# Read list of associations ----

all_assocns <- read.table(paste0(filepath, "/gwascat_hormone_associations.txt"),
                          sep = "\t", header = T, quote = "")

# For bed file, get chrX, start, end, and mapped trait
all_assocns_bed <- data.frame(chr = paste0("chr", all_assocns$CHR_ID),
                              start = all_assocns$CHR_POS - 1,
                              end = all_assocns$CHR_POS,
                              mapped_trait = all_assocns$MAPPED_TRAIT)
all_assocns_bed <- all_assocns_bed %>% distinct()

# Make sure there are no NAs
all_assocns_bed <- all_assocns_bed[complete.cases(all_assocns_bed), ]
# Make sure there are no spaces in trait names
all_assocns_bed$mapped_trait <- gsub(" ", "_", all_assocns_bed$mapped_trait)

# Make sure integers as printed as integers and not in scientific
options(scipen = 999)
write.table(all_assocns_bed,
            paste0(filepath, "/gwascat_hormone_associations_hg38.bed"), sep = "\t", 
            row.names = F, quote = F, col.names = F)
options(scipen = 0)
