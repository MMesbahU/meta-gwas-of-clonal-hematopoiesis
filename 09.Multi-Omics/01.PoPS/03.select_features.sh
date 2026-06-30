#!/bin/bash

source /broad/software/scripts/useuse

reuse Anaconda3

###### May 23, 2022
#### all meta-analysed data
# combine raw and outputs from step1

	# raw files
# while read my_gene; do cp /broad/hptmp/mesbah/ukb_chip/PoPs/output_pops/pops_step1.chr1.chr1_22.lifted_hg37.GWAMA.meta_${my_gene}.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.batch1_chr.genes.raw /broad/hptmp/mesbah/ukb_chip/PoPs/output_pops/combined_output/chr1_22.lifted_hg37.GWAMA.meta_${my_gene}.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.genes.raw; done < <(echo -e "CHIP\nDNMT3A\nTET2")

# while read my_gene; do for chr in {2..22}; do awk 'NR>2' /broad/hptmp/mesbah/ukb_chip/PoPs/output_pops/pops_step1.chr${chr}.chr1_22.lifted_hg37.GWAMA.meta_${my_gene}.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.batch${chr}_chr.genes.raw >> /broad/hptmp/mesbah/ukb_chip/PoPs/output_pops/combined_output/chr1_22.lifted_hg37.GWAMA.meta_${my_gene}.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.genes.raw; done; done < <(echo -e "CHIP\nDNMT3A\nTET2")

	# Out files
# while read my_gene; do cp /broad/hptmp/mesbah/ukb_chip/PoPs/output_pops/pops_step1.chr1.chr1_22.lifted_hg37.GWAMA.meta_${my_gene}.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.batch1_chr.genes.out /broad/hptmp/mesbah/ukb_chip/PoPs/output_pops/combined_output/chr1_22.lifted_hg37.GWAMA.meta_${my_gene}.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.genes.out; done < <(echo -e "CHIP\nDNMT3A\nTET2")

# while read my_gene; do for chr in {2..22}; do awk 'NR>1' /broad/hptmp/mesbah/ukb_chip/PoPs/output_pops/pops_step1.chr${chr}.chr1_22.lifted_hg37.GWAMA.meta_${my_gene}.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.batch${chr}_chr.genes.out >> /broad/hptmp/mesbah/ukb_chip/PoPs/output_pops/combined_output/chr1_22.lifted_hg37.GWAMA.meta_${my_gene}.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.genes.out; done; done < <(echo -e "CHIP\nDNMT3A\nTET2")

## Submit Jobs: while read genename;do qsub -R y -wd /broad/hptmp/mesbah/ukb_chip/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=60G -l h_rt=30:00:00 -N step3.POPs.${genename}.refUKB200k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/PoPs/03.select_features.sh /medpop/esp2/projects/software/PoPS/data/PoPS.features.txt.gz /broad/hptmp/mesbah/ukb_chip/PoPs/output_pops/combined_output/chr1_22.lifted_hg37.GWAMA.meta_${genename}.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies /broad/hptmp/mesbah/ukb_chip/PoPs/output_pops/combined_output/chr1_22.pops_step3.has${genename} /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/PoPs/pops/pops.feature_selection.py; done < <(echo -e "CHIP\nDNMT3A\nTET2")

##########################

###
########## 29 Dec 2021
## Combine per chr output from step2
# while read my_gene; do cp magma_out.ukb_all_200k.chr1.has${my_gene}.batch1_chr.genes.raw chr1_22.magma_out.ukb_all_200k.has${my_gene}.genes.raw; done < <(echo -e "CHIP\nDNMT3A\nTET2")

# while read my_gene; do for chr in {2..22}; do awk 'NR>2' magma_out.ukb_all_200k.chr${chr}.has${my_gene}.batch${chr}_chr.genes.raw >> chr1_22.magma_out.ukb_all_200k.has${my_gene}.genes.raw; done; done < <(echo -e "CHIP\nDNMT3A\nTET2")

## while read genename;do qsub -R y -wd /broad/hptmp/mesbah/ukb_chip/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=60G -l h_rt=30:00:00 -N step3.ukb200k_magma.${genename}.refUKB200k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/PoPs/03.select_features.sh /broad/hptmp/mesbah/ukb_chip/PoPs/pops_features/PoPS.features.txt.gz /broad/hptmp/mesbah/ukb_chip/PoPs/output_pops/chr1_22.magma_out.ukb_all_200k.has${genename} /broad/hptmp/mesbah/ukb_chip/PoPs/output_pops/chr1_22.pops_step3.has${genename} /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/PoPs/pops/pops.feature_selection.py; done < <(echo -e "CHIP\nDNMT3A\nTET2")

#####################
##
## while read genename; do qsub -R y -wd /broad/hptmp/mesbah/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=60G -l h_rt=30:00:00 -N step2.ukbEUR_magma.${genename}.ref1kg /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/PoPs/03.select_features.sh /broad/hptmp/mesbah/ukb_chip/meta_gwas/pops/data/PoPS.features.txt.gz /broad/hptmp/mesbah/ukb_chip/meta_gwas/pops/gwas/ukb/magma_out/magma_out.ukb_eur.chr1_22.has${genename} /broad/hptmp/mesbah/ukb_chip/meta_gwas/pops/gwas/ukb/magma_out/pops_step2.has${genename} /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/PoPs/pops/pops.feature_selection.py; done < <(echo -e "CHIP\nDNMT3A\nTET2")

##
pops_features=${1}

step1_output=${2}

step2_output=${3}

pops_scripts=${4} # /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/PoPs/pops/pops.feature_selection.py

# chr=${SGE_TASK_ID}

# out_prefix=${out_dir}/magma_out.ukb_eur.chr${chr}.has${genename}

######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
#########

## Run MAGMA
python ${pops_scripts} \
	--features ${pops_features} \
	--gene_results ${step1_output} \
	--out ${step2_output}

######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

