Scripts to assess extent of polygenicity and overlap between: infertility, reproductive hormones, reproductive diseases, and obesity. Only run traits with significant heritability Z>4 as polygenicity estimation is less powered than heritability estimation. Summary statistics should be formatted as LDSC format. Scripts adapted from: https://github.com/comorment/mixer/blob/main/usecases/mixer_real.md.  

**run_mixer.sh** - Run univariate and bivariate fit and test jobs. These are array jobs that run from 1-20 repeats. 
**batch_submit_mixer.R** - Submit mixer jobs for each pair of phenotypes being tested.

