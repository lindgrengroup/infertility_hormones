library(tidyverse)

gfile <- gzfile("ukb31063_h2_z4.02Oct2019.tsv.gz", 
                "rt")
df <- read.table(gfile, sep = "\t",
                 header = T, stringsAsFactors = F)

df <- df %>%
  select(all_of(c("phenotype", "description", 
                  "sex", "variable_type",
                  "n_cases", "n_controls", "h2_liability")))

# Filter to binary phenotypes with prevalence < 5% in UKB to reflect infertility prevalence

sub_df <- df %>% 
  filter(variable_type == "binary" & n_cases / (n_cases + n_controls) <= 0.05)

sum(sub_df$h2_liability < 0.1)

