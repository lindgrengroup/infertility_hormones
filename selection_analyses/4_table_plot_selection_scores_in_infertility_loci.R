# Author: Samvida S. Venkatesh
# Date: 02/10/23

library(tidyverse)
library(scales)
library(ggpubr)
library(ggrepel)
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

# Read lead SNPs ----

lead_snps <- read.table(paste0(mainpath, "/all_hormones_infertility_lead_variants.txt"),
                        sep = "\t", header = T, stringsAsFactors = F)

# Read selection scores ----

sds_scores <- read.table("/well/lindgren/samvida/Resources/selection_scores_public/SDS_UK10K_n3195_release_Sep_19_2016/hg38/SDS_UK10K_n3195_release_Sep_19_2016_with_pvals_hg38.txt",
                         sep = " ", header = T, stringsAsFactors = F)
SDS_LOW <- quantile(sds_scores$SDS, 0.025)
SDS_HIGH <- quantile(sds_scores$SDS, 0.975)

betascan_scores <- read.table("/well/lindgren/samvida/Resources/selection_scores_public/BetaScan2/hg38/all_B2std_hg38.txt",
                              sep = " ", header = T, stringsAsFactors = F)
BETA_HIGH <- quantile(betascan_scores$Beta2_std, 0.95)

# Read scores for GWASCatalog regions only ----

sds_gwascat <- read.table("/well/lindgren/samvida/hormones_infertility/selection_analyses/data/gwascat_sds_scores.txt",
                         sep = "\t", header = F, stringsAsFactors = F)
colnames(sds_gwascat) <- colnames(sds_scores)
sds_perc <- ecdf(sort(sds_gwascat$SDS))

betascan_gwascat <- read.table("/well/lindgren/samvida/hormones_infertility/selection_analyses/data/gwascat_betascan_scores.txt",
                              sep = "\t", header = F, stringsAsFactors = F)
colnames(betascan_gwascat) <- colnames(betascan_scores)
betascan_perc <- ecdf(sort(betascan_gwascat$Beta2_std))

# Function to extract scores for all variants within 10 kb of lead variant ----

WINDOW_SIZE <- 10000 # bp
getWindowScores <- function (chr_lead, pos_lead) {
  sds_window <- sds_scores %>% filter(CHR == paste0("chr", chr_lead) &
                                        POS >= pos_lead - WINDOW_SIZE &
                                        POS <= pos_lead + WINDOW_SIZE) %>%
    mutate(status = ifelse(SDS >= SDS_HIGH | SDS <= SDS_LOW, "signif", "non_sig")) %>%
    rename(BP_pos = POS, score = SDS)
  
  beta_window <- betascan_scores %>% filter(Chromosome == paste0("chr", chr_lead) &
                                              Position >= pos_lead - WINDOW_SIZE &
                                              Position <= pos_lead + WINDOW_SIZE) %>%
    mutate(status = ifelse(Beta2_std >= BETA_HIGH, "signif", "non_sig")) %>%
    rename(BP_pos = Position, score = Beta2_std)
  
  return (list(sds = sds_window,
               betascan = beta_window))
}

# Function for mini-Manhattan plots within each lead variant locus ----

col_palette <- c("#A4A4A4", "#D35C79", "#009593", "#000000")
names(col_palette) <- c("non_sig", "signif_F", "signif_M", "lead")

locusPlot <- function (locus_dat, locus_name, 
                       locus_position, 
                       scoretype = "BetaScan") {
  plot_dat <- locus_dat %>%
    mutate(status = ifelse(status == "signif" & sex_strata == "F", "signif_F",
                           ifelse(status == "signif" & sex_strata == "M", "signif_M",
                                  status)))
  # Add text for lead variant
  if (locus_position %in% plot_dat$BP_pos) {
    for_lab <- plot_dat[plot_dat$BP_pos == locus_position, ]
    for_lab$label_id <- locus_name
    for_lab$status <- "lead"
  } else {
    for_lab <- data.frame(BP_pos = locus_position,
                          score = 0, status = "lead",
                          label_id = locus_name)
  }
  plot_dat <- plot_dat %>% filter(BP_pos != locus_position)
  plot_dat <- bind_rows(plot_dat, for_lab)
  
  minpos <- min(plot_dat$BP_pos)
  maxpos <- max(plot_dat$BP_pos)
  
  if (scoretype == "BetaScan") {
    manhattan_plot <- ggplot(plot_dat, aes(x = BP_pos, y = score),
                             fill = status, colour = status) +
      geom_point(aes(fill = status, colour = status)) +
      geom_text_repel(data = for_lab, aes(label = label_id)) +
      geom_hline(yintercept = BETA_HIGH, linetype = "dashed") +
      scale_colour_manual(values = col_palette, guide = "none") +
      scale_fill_manual(values = col_palette, guide = "none") +
      scale_x_continuous(limits = c(minpos, maxpos), breaks = pretty_breaks()) +
      labs(x = "Position (hg38)", y = "Beta2_std") +
      theme(panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank()) 
  } else {
    manhattan_plot <- ggplot(plot_dat, aes(x = BP_pos, y = score),
                             fill = status, colour = status) +
      geom_point(aes(fill = status, colour = status)) +
      geom_text_repel(data = for_lab, aes(label = label_id)) +
      geom_hline(yintercept = SDS_LOW, linetype = "dashed") +
      geom_hline(yintercept = SDS_HIGH, linetype = "dashed") +
      scale_colour_manual(values = col_palette, guide = "none") +
      scale_fill_manual(values = col_palette, guide = "none") +
      scale_x_continuous(limits = c(minpos, maxpos), breaks = pretty_breaks()) +
      labs(x = "Position (hg38)", y = "SDS") +
      theme(panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank()) 
  }
  return (manhattan_plot)
}

# Apply across all hormone and infertility traits ----

hormone_selection_scores <- lapply(HORMONE_STRATA, function (hr) {
  cat(paste0("Hormone: ", hr), "\n")
  snps_to_plot <- lead_snps %>% filter(trait == hr)
  if (nrow(snps_to_plot) > 0) {
    selection_scores <- lapply(1:nrow(snps_to_plot), function (i) {
      return (getWindowScores(chr_lead = snps_to_plot$chrom[i],
                              pos_lead = snps_to_plot$pos_hg38[i]))
    })
    names(selection_scores) <- snps_to_plot$rsid
    
    # Get summary of min/max SDS and max Beta scores
    summ_dat <- lapply(snps_to_plot$rsid, function (i) { 
      # cat(paste0("i: ", i), "\n")
      sds_df <- selection_scores[[i]]$sds
      beta_df <- selection_scores[[i]]$betascan
      return (data.frame(rsid = i,
                         n_SDS = nrow(sds_df),
                         min_SDS = min(sds_df$score),
                         min_SDS_perc = sds_perc(min(sds_df$score)),
                         max_SDS = max(sds_df$score),
                         max_SDS_perc = sds_perc(max(sds_df$score)),
                         n_beta = nrow(beta_df),
                         max_beta = max(beta_df$score),
                         max_beta_perc = betascan_perc(max(beta_df$score))))
    })
    summ_dat <- bind_rows(summ_dat)
    summ_dat$trait <- hr
    
    # Get SDS plots
    SDS_plots <- lapply(snps_to_plot$rsid, function (i) {
      locus_dat <- selection_scores[[i]]$sds
      if (nrow(locus_dat) > 0) {
        locus_dat$sex_strata <- "F"
        if (grepl("_M", hr)) locus_dat$sex_strata <- "M"
        
        lplot <- locusPlot(locus_dat = locus_dat, 
                           locus_name = i, 
                           locus_position = snps_to_plot$pos_hg38[snps_to_plot$rsid == i], 
                           scoretype = "SDS")
      } else lplot <- NULL
      return (lplot)
    })
    pdf(paste0(mainpath, "/plots/SDS_scores_", hr, ".pdf"),
        onefile = T)
    print(ggarrange(plotlist = SDS_plots, 
              ncol = 2, nrow = 3, common.legend = TRUE))
    dev.off()
    
    # Get BetaScan2 plots
    beta2_plots <- lapply(snps_to_plot$rsid, function (i) {
      locus_dat <- selection_scores[[i]]$betascan
      if (nrow(locus_dat) > 0) {
        locus_dat$sex_strata <- "F"
        if (grepl("_M", hr)) locus_dat$sex_strata <- "M"
        
        lplot <- locusPlot(locus_dat = locus_dat, 
                           locus_name = i, 
                           locus_position = snps_to_plot$pos_hg38[snps_to_plot$rsid == i], 
                           scoretype = "BetaScan")
      } else lplot <- NULL
      return (lplot)
    })
    pdf(paste0(mainpath, "/plots/Beta2std_scores_", hr, ".pdf"))
    print(ggarrange(plotlist = beta2_plots, 
              ncol = 2, nrow = 3, common.legend = TRUE))
    dev.off()
  } else summ_dat <- NULL
  return (summ_dat)
})
hormone_selection_scores <- bind_rows(hormone_selection_scores)
write.table(hormone_selection_scores,
            paste0(mainpath, "/results/hormone_loci_selection_scores.txt"),
            sep = "\t", row.names = F, quote = F)

# Format for thesis table
hormone_tt <- hormone_selection_scores %>%
  mutate(write_min_sds = paste0(signif(min_SDS, 3), "(", signif(min_SDS_perc*100, 3), "%)"),
         write_max_sds = paste0(signif(max_SDS, 3), "(", signif(max_SDS_perc*100, 3), "%)"),
         write_max_beta = paste0(signif(max_beta, 3), "(", signif(max_beta_perc*100, 3), "%)")) %>%
  select(all_of(c("trait", "rsid",
                  "n_SDS", "write_min_sds", "write_max_sds",
                  "n_beta", "write_max_beta")))
write.table(hormone_tt,
            paste0(mainpath, "/results/hormone_loci_selection_scores_for_tt.txt"),
            sep = "\t", row.names = F, quote = F)


infert_selection_scores <- lapply(INFERT_STRATA, function (infert) {
  snps_to_plot <- lead_snps %>% filter(trait == infert)
  cat(paste0("Trait: ", infert), "\n")
  if (nrow(snps_to_plot) > 0) {
    selection_scores <- lapply(1:nrow(snps_to_plot), function (i) {
      return (getWindowScores(chr_lead = snps_to_plot$chrom[i],
                              pos_lead = snps_to_plot$pos_hg38[i]))
    })
    names(selection_scores) <- snps_to_plot$rsid
    
    # Get summary of min/max SDS and max Beta scores
    summ_dat <- lapply(snps_to_plot$rsid, function (i) { 
      sds_df <- selection_scores[[i]]$sds
      beta_df <- selection_scores[[i]]$betascan
      return (data.frame(rsid = i,
                         n_SDS = nrow(sds_df),
                         min_SDS = min(sds_df$score),
                         min_SDS_perc = sds_perc(min(sds_df$score)),
                         max_SDS = max(sds_df$score),
                         max_SDS_perc = sds_perc(max(sds_df$score)),
                         n_beta = nrow(beta_df),
                         max_beta = max(beta_df$score),
                         max_beta_perc = betascan_perc(max(beta_df$score))))
    })
    summ_dat <- bind_rows(summ_dat)
    summ_dat$trait <- infert
    
    # Get SDS plots
    SDS_plots <- lapply(snps_to_plot$rsid, function (i) {
      locus_dat <- selection_scores[[i]]$sds
      locus_dat$sex_strata <- "F"
      if (infert == "male_infertility") locus_dat$sex_strata <- "M"
      
      lplot <- locusPlot(locus_dat = locus_dat, 
                         locus_name = i, 
                         locus_position = snps_to_plot$pos_hg38[snps_to_plot$rsid == i], 
                         scoretype = "SDS")
      return (lplot)
    })
    pdf(paste0(mainpath, "/plots/SDS_scores_", infert, ".pdf"))
    print(ggarrange(plotlist = SDS_plots, 
              ncol = 2, nrow = 3, common.legend = TRUE))
    dev.off()
    
    # Get BetaScan2 plots
    beta2_plots <- lapply(snps_to_plot$rsid, function (i) {
      locus_dat <- selection_scores[[i]]$betascan
      locus_dat$sex_strata <- "F"
      if (infert == "male_infertility") locus_dat$sex_strata <- "M"
      
      lplot <- locusPlot(locus_dat = locus_dat, 
                         locus_name = i, 
                         locus_position = snps_to_plot$pos_hg38[snps_to_plot$rsid == i], 
                         scoretype = "BetaScan")
      return (lplot)
    })
    pdf(paste0(mainpath, "/plots/Beta2std_scores_", infert, ".pdf"))
    print(ggarrange(plotlist = beta2_plots, 
              ncol = 2, nrow = 3, common.legend = TRUE))
    dev.off()
    
  } else summ_dat <- NULL
  return (summ_dat)
})
infert_selection_scores <- bind_rows(infert_selection_scores)
write.table(infert_selection_scores,
            paste0(mainpath, "/results/infertility_loci_selection_scores.txt"),
            sep = "\t", row.names = F, quote = F)

# Format for thesis table
infert_tt <- infert_selection_scores %>%
  mutate(write_min_sds = paste0(signif(min_SDS, 3), "(", signif(min_SDS_perc*100, 3), "%)"),
         write_max_sds = paste0(signif(max_SDS, 3), "(", signif(max_SDS_perc*100, 3), "%)"),
         write_max_beta = paste0(signif(max_beta, 3), "(", signif(max_beta_perc*100, 3), "%)")) %>%
  select(all_of(c("trait", "rsid",
                  "n_SDS", "write_min_sds", "write_max_sds",
                  "n_beta", "write_max_beta")))
write.table(infert_tt,
            paste0(mainpath, "/results/infertility_loci_selection_scores_for_tt.txt"),
            sep = "\t", row.names = F, quote = F)
