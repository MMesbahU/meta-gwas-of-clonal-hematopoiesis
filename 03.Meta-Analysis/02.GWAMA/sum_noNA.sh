#!/bin/bash

source /broad/software/scripts/useuse

## while read ANC; do while read pheno; do qsub -R y -wd /broad/hptmp/mesbah/dataset/ch_gwas/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=20G -l h_rt=10:00:00 -N remove_na.${ANC}_${pheno} /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/sum_noNA.sh /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/chr1_22.MultiANC.${ANC}.has${pheno}.regenie.tsv /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb; done < <(echo -e "CH\nDNMT3A\nTET2\nCHvaf10\nASXL1\nSF\nDDR"); done < <(echo -e "ukb200k_N193342\nukb250k_N243350")

# 
inSum=${1}

outDir=${2}

#

awk 'NF == 9 && $0 !~ /NA/' ${inSum} > ${outDir}/$(basename ${inSum} ".tsv").noNA.tsv


