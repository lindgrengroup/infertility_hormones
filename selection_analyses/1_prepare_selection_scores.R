# Author: Samvida S. Venkatesh
# Date: 02/10/23

library(tidyverse)
theme_set(theme_bw())

# Convert SDS normalised scores to P-values

mainpath <- ("/well/lindgren/samvida/Resources/selection_scores_public")

getQQ <- function (obs_scores) {
  nvals <- length(obs_scores)
  # Get expected p-values
  obs_scores <- sort(obs_scores)
  exp_scores <- qnorm(p = (1:nvals - 0.5)/nvals, 
                      mean = 0, sd = 1)
  
  # CI based on the beta distribution
  # parameters: a = k, b = n-k where k=seq(1, n)
  k <- 1:nvals
  beta_lci <- qbeta(0.025, k, nvals-k)
  beta_uci <- qbeta(0.975, k, nvals-k)
  
  plot_qq <- data.frame(obs = obs_scores, exp = exp_scores,
                        lci = beta_lci, uci = beta_uci)
  
  plot_res <- ggplot(plot_qq, 
                     aes(x = exp, y = obs)) +
    geom_ribbon(aes(ymin = lci, ymax = uci), 
                fill = "grey", alpha = 0.5) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    labs(x = "Expected Z-score", y = "Observed Z-score")
  return (plot_res)
}

# Read Beta2 scores, normalise, and add p-value

beta2_scores <- lapply(CHRS, function (chr) {
  bdf <- read.table(paste0(mainpath, "/BetaScan2/StdB2Scores/chr",
                           chr, "_B2std.out"),
                    sep = "\t", header = T, stringsAsFactors = F)
  
  mean_beta2 <- mean(bdf$Beta2_std)
  sd_beta2 <- sd(bdf$Beta2_std)
  bdf$Beta2_std_zscore <- (bdf$Beta2_std - mean_beta2)/sd_beta2
  
  bdf_qq <- getQQ(bdf$Beta2_std_zscore)
  ggsave(paste0(mainpath, "/plots/BetaScan2_chr", chr, "_qq.png"), 
         bdf_qq)
  
  bdf$pval <- 2*pnorm(bdf$Beta2_std_zscore, lower.tail = F)
  write.table(bdf, paste0(mainpath, "/BetaScan2/StdB2Scores/chr",
                          chr, "_B2std_with_pvals.out"),
              sep = "\t", row.names = F, quote = F)
  return ()
})

# Read SDS scores, check that they look normal, and add p-value
sds_gz <- gzfile(paste0(mainpath, "/SDS_UK10K_n3195_release_Sep_19_2016/SDS_UK10K_n3195_release_Sep_19_2016.tab.gz"),
                 "rt")
sds_scores <- read.table(sds_gz, sep = "\t", header = T,
                         stringsAsFactors = F)
# Scores are already normal so don't worry about this
# mean_sds <- mean(sds_scores$SDS)
# sd_sds <- sd(sds_scores$SDS)
# sds_scores$SDS_zscore <- (sds_scores$SDS - mean_sds)/sd_sds

sds_qq <- getQQ(sds_scores$SDS)
ggsave(paste0(mainpath, "/plots/SDS_UK10K_n3195_qq.png"), 
       sds_qq)

sds_scores$pval <- 2*pnorm(sds_scores$SDS, lower.tail = F)
write.table(sds_scores, paste0(mainpath, "/SDS_UK10K_n3195_release_Sep_19_2016/SDS_UK10K_n3195_release_Sep_19_2016_with_pvals.txt"),
            sep = "\t", row.names = F, quote = F)


