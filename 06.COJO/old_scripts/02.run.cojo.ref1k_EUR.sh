#!/bin/bash

# for files in $(ls -l /broad/hptmp/mesbah/ukb_chip/meta_gwas/cojo/cojo_input/summary_for_cojo.lifted_hg37.GWAMA.m*.tsv | awk '{print $NF}'); do qsub -t 1-22 -R y -wd /broad/hptmp/mesbah/ukb_chip/meta_gwas/cojo/tmpdir -N run_cojo -l h_rt=20:00:00 -l h_vmem=10G -pe smp 4 -binding linear:4 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.COJO/02.run.cojo.ref1k_EUR.sh ${files} /broad/hptmp/mesbah/ukb_chip/meta_gwas/cojo/cojo_1k_EUR 0.01 0.8 5e-8 4; done

## UKB: for files in $(ls -l /broad/hptmp/mesbah/ukb_chip/meta_gwas/cojo/ukb_only/input/summary_for_cojo.chr1_22.has*.21Aug2021_ukbEUR.tsv | awk '{print $NF}'); do qsub -t 1-22 -R y -wd /broad/hptmp/mesbah/ukb_chip/meta_gwas/cojo/tmpdir -N ukb_cojo_ref1k -l h_rt=20:00:00 -l h_vmem=10G -pe smp 2 -binding linear:2 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.COJO/02.run.cojo.ref1k_EUR.sh ${files} /broad/hptmp/mesbah/ukb_chip/meta_gwas/cojo/ukb_only/output 0.01 0.8 5e-6 2; done
#### 
source /broad/software/scripts/useuse

use Bcftools
use Tabix

gcta=/medpop/esp2/mesbah/tools/gcta_1.93.2beta/gcta64

### Read command line inputs
cojo_input=${1}
cojo_output_dir=${2}

MAF=${3} # 0.01

LD_threshold=${4} # 0.8

pval_threshold=${5} # 5e-8

threads=${6} # 6

test_chr=${SGE_TASK_ID}

# LD_Reference=/broad/hptmp/mesbah/ukb10k/ukb10k_eur.imp_chr${test_chr}_v3
LD_Reference=/broad/hptmp/mesbah/ukb_chip/meta_gwas/pops/data/1000G.EUR

cojo_output=${cojo_output_dir}/ref1kg_chr${test_chr}_cojo.$(basename ${cojo_input} ".tsv")

## Run COJO
${gcta} \
	--bfile ${LD_Reference} \
	--cojo-file ${cojo_input} \
	--chr ${test_chr} \
	--maf ${MAF} \
	--cojo-collinear ${LD_threshold} \
	--cojo-p ${pval_threshold} \
	--cojo-slct \
	--cojo-wind 10000 \
	--diff-freq 0.2 \
	--thread-num ${threads} \
	--out ${cojo_output}


