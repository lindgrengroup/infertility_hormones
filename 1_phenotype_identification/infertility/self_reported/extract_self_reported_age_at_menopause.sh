#!/bin/bash 

awk -v RS='\t' '/f.3581.0.0/{print NR; exit}' /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.tab
 
awk '{print $1 "," $906 "," $907 "," $908}' /well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.tab > /well/lindgren/UKBIOBANK/laura/infertility/self_reported/self_reported_age_at_menopause_3581_0_1_2.txt

covars<-cv[,which(names(cv) %in% c("eid","X74.0.0","X21003.0.0","X54.0.0","X2003.0.0","X2724.0.0"))]

