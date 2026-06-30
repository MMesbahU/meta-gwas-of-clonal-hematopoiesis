#!/bin/bash

## qsub -R y -wd /broad/hptmp/mesbah/dataset/ch_gwas/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=20G -l h_rt=10:00:00 -N prep_data /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/prep_summary.sh

## TOPMed
# while read ANC; do while read pheno; do zcat /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv.gz > /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv; done < <(echo -e "CH\nDNMT3A\nTET2\nCHvaf10\nASXL1\nSF\nDDR"); done < <(echo -e "MultiANC\nFemale\nMale\nEUR\nAFR\nAMR")

## AoU
# while read ANC; do while read pheno; do zcat /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/gwas/aou/${ANC}/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv.gz > /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/aou/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv; done < <(echo -e "CH\nDNMT3A\nTET2\nCHvaf10\nASXL1\nSF\nDDR"); done < <(echo -e "MultiANC\nFemale\nMale\nEUR\nAFR\nAMR")

## UKBB
# while read ANC; do while read pheno; do zcat /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/step2/summary/has${pheno}.chr1_22.${ANC}.regenie.tsv.gz > /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/has${pheno}.UKB450k.chr1_22.${ANC}.regenie.tsv; done < <(echo -e "CH\nDNMT3A\nTET2\nCHvaf10\nASXL1\nSF\nDDR"); done < <(echo -e "Female\nMale\nEUR\nAFR\nAMR")

## GWAMA Summary
# UKBB
while read ANC; do while read pheno; do zcat /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/chr1_22.has${pheno}.${ANC}.ukbb450k.ukb200k_N193342.ukb250k_N243350.out.gz | awk 'BEGIN {print "SNPID\tREF\tALT\tAAF\tBETA\tSE\tP\tN\tHet_Q\tDirection"} (NR > 1 && $15 == 2) {print $1 "\t" $3 "\t" $2 "\t" $4 "\t" $5 "\t" $6 "\t" $10 "\t" $16 "\t" $13 "\t" $17}' > /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/chr1_22.has${pheno}.${ANC}.ukbb450k.ukb200k_N193342.ukb250k_N243350.tsv; done < <(echo -e "CH\nCHvaf10\nDNMT3A\nTET2\nASXL1\nSF\nDDR"); done < <(echo -e "MultiANC")

