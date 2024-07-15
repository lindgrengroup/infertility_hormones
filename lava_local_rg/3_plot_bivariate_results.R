# Author: Samvida S. Venkatesh
# Date: 15/07/2024

library(tidyverse)
theme_set(theme_bw())

mainpath <- "/well/lindgren/samvida/hormones_infertility/lava_local"

# Read data ----

fnames <- list.files(paste0(mainpath, "/bivariate_results"), pattern = "*_bivar.txt")

full_fnames <- paste0(mainpath, "/bivariate_results/", fnames)
fnames <- fnames[sapply(full_fnames, file.size) > 1] # get rid of empty files

lava_bivar <- lapply(fnames, function (fn) {
  res <- read.table(paste0(mainpath, "/bivariate_results/", fn), 
                    header = T, stringsAsFactors = F)
  return (res)
})
lava_bivar <- bind_rows(lava_bivar)

# Manhattan plot for local rG, infertility vs all other phenotypes ----

ALL_PHENOS <- unique(c(lava_bivar$phen1, lava_bivar$phen2))

INFERT <- c("female_infertility_analysis1", "female_infertility_analysis3",
            "female_infertility_analysis5")
HORMONES <- c("fsh_f", "testosterone_f", "williams_tsh", "verdiesen_amh")
REPRO <- c("endometriosis", "hmb", "pcos", "uterine_fibroids", "huber_twinning")
OBESITY <- ALL_PHENOS[!ALL_PHENOS %in% c(INFERT, HORMONES, REPRO)]

# Only retain phenotype pairs for infertility x trait 

lava_infert <- lava_bivar %>%
  filter(phen1 %in% INFERT | phen2 %in% INFERT)

# Add FDR p-value to check significance 
PTHRESH_BONF <- 0.05/nrow(lava_infert)
lava_infert$p_fdr <- p.adjust(lava_infert$p, method = "fdr")
# Closest p-value to FDR=0.05 
PTHRESH_FDR <- lava_infert$p[which.min(abs(lava_infert$p_fdr - 0.05))]

nonsig_col_palette <- c("#A9A9A9", "#D3D3D3")
names(nonsig_col_palette) <- c("odd_nonsig", "even_nonsig")

sig_col_palette <- c("#D35C79", "#FFFFCC", "#009593")
names(sig_col_palette) <- c("high", "mid", "low")

plot_dat <- lava_infert %>% 
  mutate(chrom = gsub(":.*", "", chrpos),
         chrom = as.numeric(gsub("chr", "", chrom)),
         status = ifelse(chrom %% 2 == 0 & p_fdr >= 0.05, "even_nonsig",
                         ifelse(chrom %% 2 != 0 & p_fdr >= 0.05, "odd_nonsig",
                                ifelse(p_fdr < 0.05, "signif", NA))),
         shape_class = ifelse(phen1 %in% HORMONES | phen2 %in% HORMONES, "hormone",
                              ifelse(phen1 %in% REPRO | phen2 %in% REPRO, "repro",
                                     ifelse(phen1 %in% OBESITY | phen2 %in% OBESITY, "obesity", NA)))) %>%
  # get chromosome length
  group_by(chrom) %>% summarise(chr_len = max(loc) - min(loc)) %>%
  # get chromosome position
  mutate(tot = as.numeric(cumsum(chr_len) - chr_len)) %>% select(-chr_len) %>%
  # add to original results dataset
  left_join(plot_dat, ., by = c("chrom" = "chrom")) %>%
  # add cumulative position of each locus
  arrange(chrom, loc) %>% mutate(loc_pos = as.numeric(loc + tot))

plot_dat <- plot_dat[complete.cases(plot_dat),]

# Axis should just show chromosome number
axisdf <- plot_dat %>% group_by(chrom) %>% 
  summarise(centre = (max(loc_pos) + min(loc_pos)) / 2)

# Plot
manhattan_plot <- ggplot(plot_dat, aes(x = loc_pos, y = -log10(p))) +
  geom_point(data = plot_dat %>% filter(status == "odd_nonsig"), 
             fill = "#A9A9A9", colour = "#A9A9A9",
             shape = 19, size = 1) +
  geom_point(data = plot_dat %>% filter(status == "even_nonsig"), 
             fill = "#D3D3D3", colour = "#D3D3D3",
             shape = 19, size = 1) +
  geom_point(data = plot_dat %>% filter(!status %in% c("even_nonsig", "odd_nonsig")), 
             aes(fill = rho, colour = rho,
                 shape = shape_class), size = 2) +
  # Bonferroni threshold line
  geom_hline(yintercept = -log10(PTHRESH_BONF), linetype = "dashed") +
  # FDR line
  geom_hline(yintercept = -log10(PTHRESH_FDR), linetype = "dashed") +
  scale_color_gradient2(low = "#009593",
                        mid = "#FFFFCC",
                        high = "#D35C79") +
  scale_fill_gradient2(low = "#009593",
                       mid = "#FFFFCC",
                       high = "#D35C79") +
  scale_x_continuous(label = axisdf$chrom, breaks = axisdf$centre) +
  scale_y_continuous(limits = c(0, NA)) +
  theme(panel.border = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

ggsave(paste0(mainpath, "/infert_vs_all_bivariate_manhattan.png"),
       width = 14, height = 7, units = "in", manhattan_plot)

# Choose loci for multivariate analysis follow-up ----

sig_res_follow <- lava_infert %>%
  filter(p_fdr < 0.05) 

# Write table for manuscript

sig_res_follow <- sig_res_follow %>%
  mutate(chrom = gsub(":.*", "", chrpos),
         chrom = gsub("chr", "", chrom),
         start = gsub("chr.*:", "", chrpos),
         start = gsub("-.*", "", start),
         stop = gsub(".*-", "", chrpos)) %>%
  arrange(chrom, start, phen1, phen2) %>%
  mutate(local_rg_ci = paste0(signif(rho, 3), " (", signif(rho.lower, 3), " - ", signif(rho.upper, 3)),
         local_rg_p = signif(p, 3))

sig_write <- sig_res_follow %>%
  select(all_of(c("chrpos", "phen1", "phen2", 
                  "local_rg_ci", "local_rg_p", "nsnps")))
write.table(sig_write, paste0(mainpath, "/infert_vs_all_bivariate_fdrsig_table.txt"),
            sep = "\t", row.names = F, quote = F)

# Write list of loci (chr, start, stop) for multivariate follow-up

loci_for_mult <- sig_res_follow %>%
  select(all_of(c("chrom", "start", "stop"))) %>%
  distinct()

write.table(loci_for_mult, paste0(mainpath, "/loci_for_multivariate_analysis.txt"),
            sep = "\t", row.names = F, quote = F)
