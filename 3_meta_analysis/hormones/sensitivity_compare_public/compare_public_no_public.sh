#!/bin/bash

cd /well/lindgren/samvida/hormones_infertility/meta_results_230613_no_public

cat ../with_public_strata_list.txt | while read STRATA; do
	sed -i 's/chr23:/chrX:/g' ${STRATA}_1.tbl > 
	grep -fwF compare_public/${STRATA}_lead_snps.txt ${STRATA}_1.tbl > \
	compare_public/no_public_sumstats_${STRATA}.txt
done


