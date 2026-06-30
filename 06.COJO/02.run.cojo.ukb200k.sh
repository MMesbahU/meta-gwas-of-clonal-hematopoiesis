#!/bin/bash


###

# EUR:
# while read trait; do qsub -t 1-22 -R y -wd /broad/hptmp/mesbah/gwas/cojo/tmpdir -N cojo_${trait}.eur.ref_ukb200k -l h_rt=20:00:00 -l h_vmem=20G -pe smp 5 -binding linear:5 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/06.COJO/02.run.cojo.ukb200k.sh /broad/hptmp/mesbah/gwas/cojo/input/ref_ukb200k.summary_4_cojo.lifted_hg37.eur_metaGWAS.${trait}.GWAMA.hg37_dbSNP.tsv /broad/hptmp/mesbah/gwas/cojo/ukb200k/output 0.001 0.8 5e-8 5 /broad/hptmp/mesbah/ukb_chip/ukb200k_reference/plink_ukb200k ref_ukb200k.eur_metaGWAS.${trait}; done < <(echo -e "CHIP\nDNMT3A\nTET2")

## Multi-ancestry:
## while read trait; do qsub -t 1-22 -R y -wd /broad/hptmp/mesbah/gwas/cojo/tmpdir -N cojo_${trait}.multiAncestry.ref_ukb200k -l h_rt=20:00:00 -l h_vmem=20G -pe smp 5 -binding linear:5 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/06.COJO/02.run.cojo.ukb200k.sh /broad/hptmp/mesbah/gwas/cojo/input/ref_ukb200k.summary_4_cojo.lifted_hg37.metaGWAS.${trait}.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.hg37_dbSNP.tsv /broad/hptmp/mesbah/gwas/cojo/ukb200k/output 0.001 0.8 5e-8 5 /broad/hptmp/mesbah/ukb_chip/ukb200k_reference/plink_ukb200k ref_ukb200k.multiAncestry_metaGWAS.${trait}; done < <(echo -e "CHIP\nDNMT3A\nTET2")


####### March 8, 2022
## "in_cojo_eaf001_min2Studies.8March2022.lifted_hg37.GWAMA.meta_CHIP.allsamples_topmed_ukbb_mgbb_BioVU.tsv " 
#n SNP (MAF>=0.1, n_studies>=2): 21331331
# GWAS Sig (N SNPs with P<=5e-8) =
  # FALSE     TRUE 
  # 21330185     1146 
  # FDR sig SNPs (N SNPs with FDR<=0.05)=
  #    FALSE     TRUE 
 #     21329248     2083 
#      max P with FDR<=0.05=[1] 4.87686e-06
# qsub -t 1-7 -R y -wd /broad/hptmp/mesbah/ukb_chip/tmpdir -N cojo_chip.meta4 -l h_rt=20:00:00 -l h_vmem=20G -pe smp 5 -binding linear:5 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.COJO/02.run.cojo.ukb200k.sh /broad/hptmp/mesbah/ukb_chip/cojo/cojo_input/in_cojo_eaf001_min2Studies.8March2022.lifted_hg37.GWAMA.meta_CHIP.allsamples_topmed_ukbb_mgbb_BioVU.tsv /broad/hptmp/mesbah/ukb_chip/cojo/output_cojo/allsamples 0.001 0.8 4.85233e-6 5 /broad/hptmp/mesbah/ukb_chip/ukb200k_reference/plink_ukb200k hasCHIP.cojo_eaf001_min2Studies.8March2022.lifted_hg37.GWAMA.meta_CHIP.allsamples_topmed_ukbb_mgbb_BioVU.ld8_fdr4.85233e6.Ref_ukb200k

## while read lines; do qsub -t 1-22 -R y -wd /broad/hptmp/mesbah/ukb_chip/tmpdir -N cojo_$(echo ${lines} | awk '{print $1}' ) -l h_rt=20:00:00 -l h_vmem=20G -pe smp 5 -binding linear:5 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.COJO/02.run.cojo.ukb200k.sh /broad/hptmp/mesbah/ukb_chip/cojo/cojo_input/$(echo ${lines} | awk '{print $6}' ) /broad/hptmp/mesbah/ukb_chip/cojo/output_cojo/allsamples 0.001 0.8 $(echo ${lines} | awk '{print $2}') 5 /broad/hptmp/mesbah/ukb_chip/ukb200k_reference/plink_ukb200k $(echo ${lines} | awk '{print $1}' ).cojo_eaf001_min2Studies.8March2022.ld8_fdr.Ref_ukb200k; done < <(awk 'NR>1 && NR<8' /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.COJO/)

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

#GSA_53K.merged.chr1.mgbb53k_CHIP_samples
LD_Reference=${ref_dir}/ukb_imp_chr${test_chr}_v3.ukb200k_CHIP_samples # /broad/hptmp/mesbah/ukb_chip/ukb200k_reference/plink_ukb200k/ukb_imp_chr22_v3.ukb200k_CHIP_samples

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

######
# while read traits; do head -1 /broad/hptmp/mesbah/gwas/cojo/mgbb53k/output/chr1.ref_mgb53k.eur_metaGWAS.CHIP.jma.cojo > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/cojo/chr1_22.ref_ukb200k.eur_metaGWAS.${traits}.jma.cojo && for chr in {1..22}; do awk 'NR>1' /broad/hptmp/mesbah/gwas/cojo/ukb200k/output/chr${chr}.ref_ukb200k.eur_metaGWAS.${traits}.jma.cojo >> /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/cojo/chr1_22.ref_ukb200k.eur_metaGWAS.${traits}.jma.cojo; done; done < <(echo -e "CHIP\nDNMT3A\nTET2")

# while read traits; do head -1 /broad/hptmp/mesbah/gwas/cojo/mgbb53k/output/chr1.ref_mgb53k.eur_metaGWAS.CHIP.jma.cojo > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/cojo/chr1_22.ref_ukb200k.multiAncestry_metaGWAS.${traits}.jma.cojo && for chr in {1..22}; do awk 'NR>1' /broad/hptmp/mesbah/gwas/cojo/ukb200k/output/chr${chr}.ref_ukb200k.multiAncetry_metaGWAS.${traits}.jma.cojo >> /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/cojo/chr1_22.ref_ukb200k.multiAncestry_metaGWAS.${traits}.jma.cojo; done; done < <(echo -e "CHIP\nDNMT3A\nTET2")


