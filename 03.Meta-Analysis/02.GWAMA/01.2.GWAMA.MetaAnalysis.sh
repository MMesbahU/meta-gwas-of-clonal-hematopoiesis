#!/bin/bash

################# Oct 3, 2023
## worked with: # Run config: CPU = 16; RAM = 104; Time: 5 Days

### Create GWAMA Input
# Oct 3, 2023
# ! while read pheno; do echo -e "/home/jupyter/ch_gwas/aou250k/aou_gwas_summary/chr1_22.AoU_multiAnc_2023.has${pheno}.August2023.regenie.tsv.gz\n/home/jupyter/ch_gwas/aou250k/sum4_cohort/has${pheno}.chr1_22.ukb200k.regenie.tsv.gz\n/home/jupyter/ch_gwas/aou250k/sum4_cohort/has${pheno}.chr1_22.ukb250k.regenie.tsv.gz\n/home/jupyter/ch_gwas/aou250k/sum4_cohort/topmed2021_has${pheno}.tsv.gz\n/home/jupyter/ch_gwas/aou250k/sum4_cohort/mgbb53k_gwama.chr1_22.has${pheno}.tsv.gz\n/home/jupyter/ch_gwas/aou250k/sum4_cohort/BioVU.saige_has${pheno}_results_merged_subset.tsv.gz" > /home/jupyter/ch_gwas/aou250k/gwama_meta/gwama_input.has${pheno}.AoU250k_UKBtopImp200k250k_topmed74k_mgbb53k_biovu54k.Oct2023.txt; done< <(echo -e "CH\nDNMT3A\nTET2")
# success!

## Run GWAMA Meta-Analysis
## Run GWAMA meta analysis in the terminal
### Run GWAMA: fixed-effect IVW
# while read pheno; do bash /home/jupyter/ch_gwas/aou250k/meta5/01.2.GWAMA.MetaAnalysis.sh /home/jupyter/ch_gwas/aou250k/gwama_meta/gwama_input.has${pheno}.AoU250k_UKBtopImp200k250k_topmed74k_mgbb53k_biovu54k.Oct2023.txt /home/jupyter/ch_gwas/aou250k/gwama_meta/metagwas.has${pheno}.AoU250k_UKBtopImp200k250k_topmed74k_mgbb53k_biovu54k.Oct2023 1>/home/jupyter/ch_gwas/aou250k/gwama_meta/log.gwamameta.has${pheno}.log 2>/home/jupyter/ch_gwas/aou250k/gwama_meta/err.gwamameta.has${pheno}.err; done < <(echo -e "CH\nDNMT3A\nTET2")
# 



################

#############
## bash 
## while read pheno; do bash /home/jupyter/ch_gwas/aou250k/meta5/01.2.GWAMA.MetaAnalysis.sh /home/jupyter/ch_gwas/aou250k/meta5/gwama_input.has${pheno}.AoU250k_topmed74k_topImpukb200k250k_mgbb53k_biovu54k.Aug2023.txt /home/jupyter/ch_gwas/aou250k/meta5/metagwas.has${pheno}.AoU250k_topmed74k_topImpukb200k250k_mgbb53k_biovu54k.Aug2023 1>/home/jupyter/ch_gwas/aou250k/meta5/gwama.has${pheno}.log 2>/home/jupyter/ch_gwas/aou250k/meta5/gwama.has${pheno}.err; done < <(echo -e "CH\nDNMT3A\nTET2")

## ### Create GWAMA Input

# while read pheno; do echo -e "/home/jupyter/ch_gwas/aou250k/aou_gwas_summary/chr1_22.AoU_multiAnc_2023.has${pheno}.August2023.regenie.tsv\n/home/jupyter/ch_gwas/aou250k/meta4/metaGWAS.has${pheno}.topmed74k_topImpukb200k250k_mgbb53k_biovu54k.Aug2023.eaf001.tsv" > /home/jupyter/ch_gwas/aou250k/meta5/gwama_input.has${pheno}.AoU250k_topmed74k_topImpukb200k250k_mgbb53k_biovu54k.Aug2023.txt; done< <(echo -e "CH\nDNMT3A\nTET2")

## Get GWAS summary: MAF>0.10% and P<5e-5
## while read pheno; do zcat metagwas.has${pheno}.AoU250k_topmed74k_topImpukb200k250k_mgbb53k_biovu54k.Aug2023.out.gz | awk '(NR==1){print $0}(NR>1 && $10<5e-5 && $4>=0.001 && $4<=0.999){print $0}' | gzip -c > p5e5_maf001.metagwas.has${pheno}.AoU250k_topmed74k_topImpukb200k250k_mgbb53k_biovu54k.Aug2023.tsv.gz; done < <(echo -e "CH\nDNMT3A\nTET2") &

######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)

#########


##
gwas_list=${1}

output_prefix=${2}

##
GWAMA=/home/jupyter/workspaces/detectionofclonalhematopoiesisofindeterminatepotentialchip/GWAMA_v2.2.2/GWAMA
##
${GWAMA} \
	-i ${gwas_list} \
	-qt \
	--name_marker SNPID \
	--name_n N \
	--name_ea ALT \
	--name_nea REF \
	--name_eaf AAF \
	--name_beta BETA \
	--name_se SE \
	--indel_alleles \
	-o ${output_prefix}

## Compress
gzip -f ${output_prefix}.out

######### Clock time #########
echo -e "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'

######### Clock time ########


