# Author: Samvida S. Venkatesh
# Date: 02/10/23

library(tidyverse)
library(scales)
theme_set(theme_bw())

mainpath <- "/well/lindgren/samvida/hormones_infertility/selection_analyses"
HORMONE_STRATA <- c("FSH_F", "FSH_M",
                    "LH_F", "LH_M",
                    "Oestradiol_F", "Oestradiol_M",
                    "Testosterone_F", "Testosterone_M")

INFERT_STRATA <- c("female_infertility_analysis1",
                   "female_infertility_analysis2",
                   "female_infertility_analysis3",
                   "female_infertility_analysis4",
                   "female_infertility_analysis5",
                   "male_infertility")

# Read selection loci ----

selection_loci <- read.table(paste0(mainpath, "/selection_loci.txt"),
                             header = T, stringsAsFactors = F)
PTHRESH <- 0.05/nrow(selection_loci)

# Function for mini-Manhattan plots within each selection locus ----

col_palette <- c("#A4A4A4", "#D35C79", "#009593")
names(col_palette) <- c("non_sig", "signif_F", "signif_M")

locusPlot <- function (locus_dat, locus_name) {
  plot_dat <- locus_dat %>%
    mutate(status = ifelse(PVAL <= PTHRESH & sex_strata == "F", "signif_F",
                           ifelse(PVAL <= PTHRESH & sex_strata == "M", "signif_M",
                                  ifelse(PVAL > PTHRESH, "non_sig", NA))))
  minpos <- min(plot_dat$BP_pos)
  maxpos <- max(plot_dat$BP_pos)
  
  manhattan_plot <- ggplot(plot_dat, aes(x = BP_pos, y = -log10(PVAL)),
                           fill = status, colour = status) +
    facet_wrap(~trait, ncol = 2, scales = "free") +
    geom_point(aes(fill = status, colour = status)) +
    geom_hline(yintercept = -log10(PTHRESH), linetype = "dashed") +
    scale_colour_manual(values = col_palette, guide = "none") +
    scale_fill_manual(values = col_palette, guide = "none") +
    scale_x_continuous(limits = c(minpos, maxpos), breaks = pretty_breaks()) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(x = "Position (hg38)", y = "-log10(Pval)") +
    theme(panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank()) +
    ggtitle(locus_name)
  return (manhattan_plot)
}

PTHRESH <- 1e-6

locusPlotforThesis <- function (locus_dat) {
  plot_dat <- locus_dat %>%
    mutate(status = ifelse(PVAL <= PTHRESH & sex_strata == "F", "signif_F",
                           ifelse(PVAL <= PTHRESH & sex_strata == "M", "signif_M",
                                  ifelse(PVAL > PTHRESH, "non_sig", NA))))
  minpos <- min(plot_dat$BP_pos)
  maxpos <- max(plot_dat$BP_pos)
  
  manhattan_plot <- ggplot(plot_dat, aes(x = BP_pos, y = -log10(PVAL)),
                           fill = status, colour = status) +
    geom_point(aes(fill = status, colour = status)) +
    geom_hline(yintercept = -log10(PTHRESH), linetype = "dashed") +
    scale_colour_manual(values = col_palette, guide = "none") +
    scale_fill_manual(values = col_palette, guide = "none") +
    scale_x_continuous(limits = c(minpos, maxpos), breaks = pretty_breaks(n = 3)) +
    scale_y_continuous(limits = c(0, NA)) +
    theme(panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(size = 6, family = "Arial")) 
  return (manhattan_plot)
}

# Apply across all loci and hormones ----

pdf(paste0(mainpath, "/plots/selection_loci_hormone_pvals.pdf"))
hormone_locus_dat <- lapply(1:nrow(selection_loci), function (i) {
  locus_name <- paste0("chr", selection_loci$CHR[i], 
                       "_", selection_loci$START_hg38[i], "_", 
                       selection_loci$END_hg38[i])
  
  hormone_dat <- lapply(HORMONE_STRATA, function (hr) {
    locus_res <- read.table(paste0(mainpath, "/data/", 
                                   locus_name, "_",
                                   hr, "_EUR.txt"),
                            sep = "\t", header = T, stringsAsFactors = F)
    
    locus_res <- locus_res[, c("ID", "GENPOS", "PVALUE")]
    colnames(locus_res) <- c("ID", "BP_pos", "PVAL")
    locus_res$trait <- hr
    locus_res$sex_strata <- "F"
    if (grepl("_M", hr)) locus_res$sex_strata <- "M"
    return (locus_res)
  })
  hormone_dat <- bind_rows(hormone_dat)
  print(locusPlot(hormone_dat, locus_name))
  
  # Create table for results with minimum p-value for each locus per strata
  summ_dat <- hormone_dat %>%
    group_by(trait) %>%
    summarise(nsnps = n(),
              min_pval = min(PVAL))
  summ_dat$locus_name <- locus_name
  return (summ_dat)
})
dev.off()

hormone_locus_dat <- bind_rows(hormone_locus_dat)
hormone_locus_wide <- hormone_locus_dat %>%
  pivot_wider(id_cols = locus_name,
              names_from = trait,
              values_from = c(nsnps, min_pval))
write.table(hormone_locus_wide, 
            paste0(mainpath, "/results/selection_loci_min_hormone_pvals.txt"),
            sep = "\t", row.names = F, quote = F)

# Apply across all loci and infertility ----

pdf(paste0(mainpath, "/plots/selection_loci_infertility_pvals.pdf"))
infert_locus_dat <- lapply(1:nrow(selection_loci), function (i) {
  locus_name <- paste0("chr", selection_loci$CHR[i], 
                       "_", selection_loci$START_hg38[i], "_", 
                       selection_loci$END_hg38[i])
  
  infert_dat <- lapply(INFERT_STRATA, function (inf) {
    locus_res <- read.table(paste0(mainpath, "/data/", 
                                   locus_name, "_",
                                   inf, "_EUR.txt"),
                            sep = "\t", header = T, stringsAsFactors = F)
    
    locus_res <- locus_res[, c("MarkerName", "pos", "P.value")]
    colnames(locus_res) <- c("ID", "BP_pos", "PVAL")
    locus_res$trait <- inf
    locus_res$sex_strata <- "F"
    if (inf == "male_infertility") locus_res$sex_strata <- "M"
    return (locus_res)
  })
  infert_dat <- bind_rows(infert_dat)
  print(locusPlot(infert_dat, locus_name))
  
  # Create table for results with minimum p-value for each locus per strata
  summ_dat <- infert_dat %>%
    group_by(trait) %>%
    summarise(nsnps = n(),
              min_pval = min(PVAL))
  summ_dat$locus_name <- locus_name
  return (summ_dat)
})
dev.off()

infert_locus_dat <- bind_rows(infert_locus_dat)
infert_locus_wide <- infert_locus_dat %>%
  pivot_wider(id_cols = locus_name,
              names_from = trait,
              values_from = c(nsnps, min_pval))
write.table(infert_locus_wide, 
            paste0(mainpath, "/results/selection_loci_min_infertility_pvals.txt"),
            sep = "\t", row.names = F, quote = F)

# For thesis (aDNA testosterone only) ----

selection_loci <- selection_loci %>% filter(Method == "aDNA")
HORMONE_STRATA <- c("Testosterone_F", "Testosterone_M")

hormone_locus_dat <- lapply(1:nrow(selection_loci), function (i) {
  locus_name <- paste0("chr", selection_loci$CHR[i], 
                       "_", selection_loci$START_hg38[i], "_", 
                       selection_loci$END_hg38[i])
  
  hormone_dat <- lapply(HORMONE_STRATA, function (hr) {
    locus_res <- read.table(paste0(mainpath, "/data/", 
                                   locus_name, "_",
                                   hr, "_EUR.txt"),
                            sep = "\t", header = T, stringsAsFactors = F)
    
    locus_res <- locus_res[, c("ID", "GENPOS", "PVALUE")]
    colnames(locus_res) <- c("ID", "BP_pos", "PVAL")
    locus_res$trait <- hr
    locus_res$sex_strata <- "F"
    if (grepl("_M", hr)) locus_res$sex_strata <- "M"
    
    ggsave(paste0(mainpath, "/plots/adna_", locus_name, "_", 
                  hr, ".png"),
           locusPlotforThesis(locus_res),
           units = "cm", height = 3.5, width = 5.6)
    return ()
  })
  return()
})

