Scripts to assess extent of polygenicity and overlap between: infertility, reproductive hormones, reproductive diseases, and obesity. Only run traits with significant heritability Z>4 as polygenicity estimation is less powered than heritability estimation. Summary statistics should be formatted as LDSC format. Scripts adapted from: https://github.com/comorment/mixer/blob/main/usecases/mixer_real.md.  

**run_mixer_fit1.sh** - Run univariate fit and test jobs for all phenotypes. These are array jobs that run from 1-20 repeats, runtime approx. 1 hour per repeat per phenotype. 
**run_mixer_fit2.sh** - Run bivariate fit and test jobs for all pairs of phenotypes, with infertility as the target of interest across different phenotypes. These are array jobs that run from 1-20 repeats, runtime approx. 2.5 hours per repeat per pair.

**batch_submit_mixer.R** - Submit fit1 and fit2 mixer jobs.

