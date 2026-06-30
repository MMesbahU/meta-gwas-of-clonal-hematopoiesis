#!/bin/bash

source /broad/software/scripts/useuse

# use Tabix

##### 
###################### MultiANC Meta-GWAS ####################
## MultiANC v1: 
# CH | DNMT3A | TET2 | CHvaf10 | ASXL1 | SF | DDR : UKB450k + AoU + TOPMed + MGBB
# while read ANC; do while read pheno; do echo -e "/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/chr1_22.has${pheno}.${ANC}.ukbb450k.ukb200k_N193342.ukb250k_N243350.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/aou/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mgbb/has${pheno}.chr1_22.${ANC}.regenie.tsv" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb450k_AoU250k_TOPMed72k_MGBB53k.txt; qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/tmpdir -pe smp 4 -binding linear:4 -l h_vmem=30G -l h_rt=10:00:00 -N ${ANC}_${pheno}.ukb_aou_top_mgb /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/01.GWAMA.Meta_Analysis.all.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb450k_AoU250k_TOPMed72k_MGBB53k.txt /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/chr1_22.has${pheno}.${ANC}.UKBB450k_AoU250k_TOPMed72k_MGBB53k; done < <(echo -e "CH\nCHvaf10\nDNMT3A\nTET2\nASXL1\nSF\nDDR"); done < <(echo -e "MultiANC")

# CH | DNMT3A | TET2: UKB450k + AoU + TOPMed + MGBB + BioVU
# while read ANC; do while read pheno; do echo -e "/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/chr1_22.has${pheno}.${ANC}.ukbb450k.ukb200k_N193342.ukb250k_N243350.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/aou/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mgbb/has${pheno}.chr1_22.${ANC}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/biovu/BioVU.saige_has${pheno}_results_merged_subset.tsv" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb450k_AoU250k_TOPMed72k_MGBB53k_BioVU54k.txt; qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/tmpdir -pe smp 4 -binding linear:4 -l h_vmem=30G -l h_rt=10:00:00 -N ${ANC}_${pheno}.ukb_aou_top_mgb_biovu /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/01.GWAMA.Meta_Analysis.all.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb450k_AoU250k_TOPMed72k_MGBB53k_BioVU54k.txt /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/chr1_22.has${pheno}.${ANC}.UKBB450k_AoU250k_TOPMed72k_MGBB53k_BioVU54k; done < <(echo -e "CH\nDNMT3A\nTET2"); done < <(echo -e "MultiANC")

## MultiANC v2: 
# CH | DNMT3A | TET2: UKB200k + UKB250k + AoU + TOPMed + MGBB + BioVU 

# CHvaf10 | ASXL1 | SF | DDR: UKB200k + UKB250k + AoU + TOPMed + MGBB

## MultiANC v3: No UKB
# CH | DNMT3A | TET2: AoU + TOPMed + MGBB + BioVU 
# CHvaf10 | ASXL1 | SF | DDR: AoU + TOPMed + MGBB

#################################################################
## Input for GWAMA
# while read pheno; do while read ANC; do echo -e "/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/has${pheno}.UKB450k.chr1_22.${ANC}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/aou/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mgbb/has${pheno}.chr1_22.${ANC}.regenie.tsv" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb450k_AoU250k_TOPMed72k_MGBB53k.txt; done < <(echo -e "AFR\nAMR\nEUR\nFemale\nMale"); done < <(echo -e "CH\nCHvaf10\nDNMT3A\nTET2\nASXL1\nSF\nDDR")

## EUR with bioVU
# while read pheno; do while read ANC; do echo -e "/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/has${pheno}.UKB450k.chr1_22.${ANC}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/aou/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mgbb/has${pheno}.chr1_22.${ANC}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/biovu/BioVU.saige_has${pheno}_results_merged_subset.tsv" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb450k_AoU250k_TOPMed72k_MGBB53k_BioVU54k.txt; done < <(echo -e "EUR"); done < <(echo -e "CH\nDNMT3A\nTET2")
## MultiANC UKB 200k 250k
# while read ANC; do while read pheno; do echo -e "/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/chr1_22.MultiANC.ukb200k_N193342.has${pheno}.regenie.noNA.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/chr1_22.MultiANC.ukb250k_N243350.has${pheno}.regenie.noNA.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/aou/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mgbb/has${pheno}.chr1_22.${ANC}.regenie.tsv" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb200k_ukb250k_AoU250k_TOPMed72k_MGBB53k.txt; done < <(echo -e "CH\nCHvaf10\nDNMT3A\nTET2\nASXL1\nSF\nDDR"); done < <(echo -e "MultiANC")
#  MultiANC UKB 200k 250k with BioVU
# while read ANC; do while read pheno; do echo -e "/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/chr1_22.MultiANC.ukb200k_N193342.has${pheno}.regenie.noNA.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/chr1_22.MultiANC.ukb250k_N243350.has${pheno}.regenie.noNA.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/aou/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mgbb/has${pheno}.chr1_22.${ANC}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/biovu/BioVU.saige_has${pheno}_results_merged_subset.tsv" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.ukbb200k_ukb250k_AoU250k_TOPMed72k_MGBB53k_BioVU54k.txt; done < <(echo -e "CH\nDNMT3A\nTET2"); done < <(echo -e "MultiANC")
## No UKBB
# while read pheno; do while read ANC; do echo -e "/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/aou/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mgbb/has${pheno}.chr1_22.${ANC}.regenie.tsv" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.AoU250k_TOPMed72k_MGBB53k.noukbb.txt; done < <(echo -e "AFR\nAMR\nEUR\nFemale\nMale\nMultiANC"); done < <(echo -e "CH\nCHvaf10\nDNMT3A\nTET2\nASXL1\nSF\nDDR")

## EUR with bioVU
# while read pheno; do while read ANC; do echo -e "/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/aou/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mgbb/has${pheno}.chr1_22.${ANC}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/biovu/BioVU.saige_has${pheno}_results_merged_subset.tsv" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.AoU250k_TOPMed72k_MGBB53k_BioVU54k.noukbb.txt; done < <(echo -e "EUR"); done < <(echo -e "CH\nDNMT3A\nTET2")

# while read pheno; do while read ANC; do echo -e "/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/aou/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mgbb/has${pheno}.chr1_22.${ANC}.regenie.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/biovu/BioVU.saige_has${pheno}_results_merged_subset.tsv" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.${ANC}.AoU250k_TOPMed72k_MGBB53k_BioVU54k.noukbb.txt; done < <(echo -e "MultiANC"); done < <(echo -e "CH\nDNMT3A\nTET2")

# Meta-Ancestry: EUR + AFR 

## EUR + AFR + AMR

## AFR + AMR

## Male + Female


###################### Female/Male Meta-GWAS ####################
## Female v1: 
# CH | DNMT3A | TET2 | CHvaf10 | ASXL1 | SF: UKB450k + AoU + TOPMed + MGBB

## Male v1: 
# CH | DNMT3A | TET2 | CHvaf10 | ASXL1 | SF: UKB450k + AoU + TOPMed + MGBB
###################################################################


###################### AFR Meta-GWAS ####################
## AFR v1: 
# CH | DNMT3A | TET2 | CHvaf10 | ASXL1 | SF: UKB450k + AoU + TOPMed + MGBB

## AFR v2: no ukb
# CH | DNMT3A | TET2 | CHvaf10 | ASXL1 | SF: AoU + TOPMed + MGBB
###################################################################


###################### AMR Meta-GWAS ####################
## AMR v1: 
# CH | DNMT3A | TET2 | CHvaf10 | ASXL1 | SF: UKB450k + AoU + TOPMed + MGBB

## AMR v2: no ukb
# CH | DNMT3A | TET2 | CHvaf10 | ASXL1 | SF: AoU + TOPMed + MGBB
###################################################################

###################### EUR Meta-GWAS ####################
## EUR v1: 
# CH | DNMT3A | TET2 | CHvaf10 | ASXL1 | SF: UKB450k + AoU + TOPMed + MGBB

# CH | DNMT3A | TET2 : UKB450k + AoU + TOPMed + MGBB + BioVU

## EUR v2: no ukb
# CH | DNMT3A | TET2 | CHvaf10 | ASXL1 | SF: AoU + TOPMed + MGBB
###################################################################


#################################################################

## UKB200k + UKB250k: UKB450k
# while read pheno; do echo -e "/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/chr1_22.MultiANC.ukb200k_N193342.has${pheno}.regenie.noNA.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/chr1_22.MultiANC.ukb250k_N243350.has${pheno}.regenie.noNA.tsv" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.MultiANC.ukbb450k.ukb200k_N193342.ukb250k_N243350.txt; qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/tmpdir -pe smp 4 -binding linear:4 -l h_vmem=30G -l h_rt=10:00:00 -N MultiANC.ukb200k_250k.${pheno} /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/01.GWAMA.Meta_Analysis.ukb450k.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.MultiANC.ukbb450k.ukb200k_N193342.ukb250k_N243350.txt /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/chr1_22.has${pheno}.MultiANC.ukbb450k.ukb200k_N193342.ukb250k_N243350; done < <(echo -e "CH\nCHvaf10\nDNMT3A\nTET2\nASXL1\nSF\nDDR")

############### Task: Worked
## 1. Run Meta-Analysis
###########################

######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
##################################
# zcat /medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/chr1_22.has${traits}.21Aug2021_ukbEUR.tsv.gz | awk '(NR==1){print "SNPID\tN\tREF\tALT\tAAF\tBETA\tSE\tP"}(NR>1 && $2~/^[0-9]+/){print $1"\t"($14+$18)"\t"$4"\t"$5"\t"$13"\t"$22"\t"$23"\t"$12}' | awk '!seen[$1]++'| gzip -c > /medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/chr1_22.ukb200k.hg37_eur_${traits}.tsv.gz
##
## 
gwas_list=${1}

output_prefix=${2}

# outDir=${3} # summary tsv file for plotting

## Run GWAMA
GWAMA=/medpop/esp2/mesbah/tools/GWAMA_v2.2.2/GWAMA 
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

## compress
gzip ${output_prefix}.out


###########


######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

