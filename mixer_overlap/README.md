Scripts to assess extent of polygenicity and overlap between: infertility, reproductive hormones, reproductive diseases, and obesity. Only run traits with significant heritability Z>4 as polygenicity estimation is less powered than heritability estimation. Summary statistics should be formatted as LDSC format. Scripts adapted from: https://github.com/comorment/mixer/blob/main/usecases/mixer_real.md.  

**run_mixer_fit1.sh** - Run univariate fit and test jobs for all phenotypes. These are array jobs that run from 1-20 repeats. 
**run_mixer_fit2.sh** - Run bivariate fit and test jobs for all pairs of phenotypes, with infertility as the target of interest (3) across different phenotypes. These are array jobs that run from 1-20 repeats.

**batch_submit_mixer.R** - Submit fit1 and fit2 mixer jobs.

