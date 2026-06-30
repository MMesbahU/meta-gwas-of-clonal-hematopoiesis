#!/bin/bash


### EUR:
# while read trait; do qsub -t 1-22 -R y -wd /broad/hptmp/mesbah/gwas/cojo/tmpdir -N cojo_${trait}.eur.ref_mgbb53k -l h_rt=20:00:00 -l h_vmem=20G -pe smp 5 -binding linear:5 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/06.COJO/02.run.cojo.mgbb53k.sh /broad/hptmp/mesbah/gwas/cojo/input/mgbb53k_ref.summary_4_cojo.eur_metaGWAS.${trait}.GWAMA.hg38_dbSNP.tsv /broad/hptmp/mesbah/gwas/cojo/mgbb53k/output 0.001 0.8 5e-8 5 /broad/hptmp/mesbah/RefSeq/mgbb53k ref_mgb53k.eur_metaGWAS.${trait}; done < <(echo -e "CHIP\nDNMT3A\nTET2")

# multi-ancestry
# while read trait; do qsub -t 1-22 -R y -wd /broad/hptmp/mesbah/gwas/cojo/tmpdir -N cojo_${trait}.multiAncestry.ref_mgbb53k -l h_rt=20:00:00 -l h_vmem=20G -pe smp 5 -binding linear:5 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/06.COJO/02.run.cojo.mgbb53k.sh /broad/hptmp/mesbah/gwas/cojo/input/mgbb53k_ref.summary_4_cojo.metaGWAS.${trait}.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.hg38_dbSNP.tsv /broad/hptmp/mesbah/gwas/cojo/mgbb53k/output 0.001 0.8 5e-8 5 /broad/hptmp/mesbah/RefSeq/mgbb53k ref_mgb53k.multiAncestry_metaGWAS.${trait}; done < <(echo -e "CHIP\nDNMT3A\nTET2")

###############
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

ref_dir=${7}

outPrefix=${8}

test_chr=${SGE_TASK_ID}

#/broad/hptmp/mesbah/RefSeq/mgbb53k/GSA_53K.merged.chr1.mgbb53k_CHIP_samples
LD_Reference=${ref_dir}/GSA_53K.merged.chr${test_chr}.mgbb53k_CHIP_samples # /broad/hptmp/mesbah/ukb_chip/ukb200k_reference/plink_ukb200k/ukb_imp_chr22_v3.ukb200k_CHIP_samples

cojo_output=${cojo_output_dir}/chr${test_chr}.${outPrefix}

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

####
## while read traits; do head -1 /broad/hptmp/mesbah/gwas/cojo/mgbb53k/output/chr1.ref_mgb53k.multiAncestry_metaGWAS.CHIP.jma.cojo > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/cojo/chr1_22.ref_mgb53k.multiAncestry_metaGWAS.${traits}.jma.cojo && for chr in {1..22}; do awk 'NR>1' /broad/hptmp/mesbah/gwas/cojo/mgbb53k/output/chr${chr}.ref_mgb53k.multiAncestry_metaGWAS.${traits}.jma.cojo >> /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/cojo/chr1_22.ref_mgb53k.multiAncestry_metaGWAS.${traits}.jma.cojo; done; done < <(echo -e "CHIP\nDNMT3A\nTET2") &
## while read traits; do head -1 /broad/hptmp/mesbah/gwas/cojo/mgbb53k/output/chr1.ref_mgb53k.eur_metaGWAS.CHIP.jma.cojo > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/cojo/chr1_22.ref_mgb53k.eur_metaGWAS.${traits}.jma.cojo && for chr in {1..22}; do awk 'NR>1' /broad/hptmp/mesbah/gwas/cojo/mgbb53k/output/chr${chr}.ref_mgb53k.eur_metaGWAS.${traits}.jma.cojo >> /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/cojo/chr1_22.ref_mgb53k.eur_metaGWAS.${traits}.jma.cojo; done; done < <(echo -e "CHIP\nDNMT3A\nTET2") &

