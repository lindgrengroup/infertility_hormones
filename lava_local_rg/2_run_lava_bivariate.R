# Author: Samvida S. Venkatesh
# Date: 08/07/2024
# Adapted from: https://github.com/josefin-werme/LAVA

library(argparse)
library(tidyverse)
library(LAVA)

parser <- ArgumentParser()
parser$add_argument("--pheno1", required=TRUE,
                    help = "First phenotype in pair")
parser$add_argument("--pheno2", required=TRUE,
                    help = "Second phenotype in pair")
parser$add_argument("--chrom", required = TRUE,
                    help = "Chromosome for loci block")
args <- parser$parse_args()

mainpath <- "/well/lindgren/samvida/hormones_infertility/lava_local"

# Read data ----

# summary statistics
input <- process.input(input.info.file = paste0(mainpath, "/input_sumstats.txt"),
                       sample.overlap.file = paste0(mainpath, "/cross_trait_sample_overlap_from_ldsc.txt"),
                       ref.prefix = "/well/lindgren/samvida/Resources/1000Genomes/for_lava/g1000_eur", # 1000 genomes on GrCh37
                       phenos = c(args$pheno1, args$pheno2)) 

# locus file
loci <- read.loci(paste0(mainpath, "/locus_files/LAVA_s2500_m25_f1_w200_chr", 
                         args$chrom, ".loci"))
NLOCI <- nrow(loci) 

# Filter out loci with no univariate genetic signals for this pair of traits ----

# Since we test N loci, P < 0.05/2495 for univariate test
PTHRESH <- 0.05/2495

bivar_res_all <- lapply(1:NLOCI, function (li) {
  res_row <- NULL
  # Some loci have negative variance estimates, which we want to ignore
  tryCatch(
    {
      locus <- process.locus(loci[li,], input)
      # only run bivariate test for loci that pass the univariate threshold
      bivar_res <- run.univ.bivar(locus, univ.thresh = PTHRESH)
      if (!is.null(bivar_res$bivar)) {
        res_row <- as.data.frame(bivar_res$bivar)
        res_row$loc <- li
        res_row$chrpos <- paste0("chr", locus$chr, ":", locus$start, "-", locus$stop)
        res_row$nsnps <- locus$n.snps
      } 
    },
    error = function(e) {
      return ( NULL )
    }
  )
  return (res_row)
})
bivar_res_all <- bind_rows(bivar_res_all)

write.table(bivar_res_all, paste0(mainpath, "/bivariate_results/",
                                  args$pheno1, "_", args$pheno2, "_chr", args$chrom,
                                  "_bivar.txt"),
            sep = "\t", row.names = F, quote = F)
