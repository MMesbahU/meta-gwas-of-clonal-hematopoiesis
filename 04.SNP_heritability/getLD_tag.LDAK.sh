#!/bin/bash

## wget https://github.com/dougspeed/LDAK/raw/refs/heads/main/ldak6.linux

mkdir -p /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/{ldak_thin,bld_ldak}

## 
## for SumHer: BLD_LDAK
## For each heritability model, there are four versions: GBR (computed using 2000 white British individuals), SAS (4214 Indian and Pakistani individuals), EAS (1279 Chinese individuals) and AFR (2577 African individuals). Click here to see a principal component plot illustrating the four different populations.

# For each population, we provide sets of pre-computed files corresponding to two SNP subsets: 1.0-1.2M non-ambiguous HapMap3 SNPs and 320-580k non-ambiguous directly genotyped SNPs. You should use whichever version best matches your summary statistics (in general, you should use the HapMap3 version if you have summary statistics from a GWAS that used imputation, otherwise use the directly genotyped version).

wget -P /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/bld_ldak/ https://genetics.ghpc.au.dk/doug/bld.ldak.{hapmap,genotyped}.{gbr,sas,eas,afr}.tagging.gz

gzip -d /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/bld_ldak/bld.ldak.{hapmap,genotyped}.{gbr,sas,eas,afr}.tagging.gz

## for genetic cor
# wget -P /broad/hptmp/mesbah/gwas/sumher/ref/ldak_thin/ https://genetics.ghpc.au.dk/doug/ldak.thin.genotyped.{gbr,sas,eas,afr}.tagging.gz

# gzip -d /broad/hptmp/mesbah/gwas/sumher/ref/ldak_thin/ldak.thin.genotyped.{gbr,sas,eas,afr}.tagging.gz

##
wget -P /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref https://genetics.ghpc.au.dk/doug/bld.zip

unzip /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/bld.zip bldnames

rm /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/bld.zip

echo -e "65 LDAK_Weightings\n66 Base_Category" >> bldnames 

