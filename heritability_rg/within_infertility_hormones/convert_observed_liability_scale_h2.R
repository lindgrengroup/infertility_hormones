# Author: Samvida S. Venkatesh
# Date: 13/11/23

#K=pop prevalence
#P=proportion of cases in study
#hsq=Heritability estimate (on observed scale)
#bigT = liability threshold
#tau = density of gaussian

INFERT <- c("female_infertility_analysis1",
            "female_infertility_analysis2",
            "female_infertility_analysis3",
            "female_infertility_analysis4",
            "female_infertility_analysis5",
            "male_infertility")

K = c(0.073106961,
      0.02523914,
      0.011282383,
      0.039290659,
      0.024979308,
      0.017177244)
names(K) <- INFERT
P = c(0.073106961,
      0.02523914,
      0.011282383,
      0.039290659,
      0.024979308,
      0.017177244)
names(P) <- INFERT
h2 = c(0.0066,
       0.0022,
       7.50E-03,
       2.20E-03,
       0.0048,
       0.0012)
names(h2) <- INFERT

h2_liab <- lapply(INFERT, function (i) {
  zv <- dnorm(qnorm(K[i]))
  
  res <- h2[i] * K[i]^2 * ( 1 - K[i])^2 / P[i] / (1-P[i]) / zv^2
  return (res)
})
h2_liab <- unlist(h2_liab)
names(h2_liab) <- INFERT

