#!/bin/bash

source /broad/software/scripts/useuse

use Tabix

############### Task: Worked
## 1. Run Meta-Analysis
## 2. sorted summary stats MAF 0.1 and N_studies>1
###########################

######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
#########

########################################## Meta4: TOPMed74k + UKB450k + MGB53k + BioVU54k ###############
# while read pheno; do echo -e "/medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2021/topmed2021_has${pheno}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/inputUKB/2023_ukbb/rap_gwas/has${pheno}.chr1_22.ukb200k.regenie.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/inputUKB/2023_ukbb/rap_gwas/has${pheno}.chr1_22.ukb250k.regenie.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/MGBB53k/mgbb53k_gwama.chr1_22.has${pheno}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_has${pheno}_results_merged_subset.tsv.gz" > /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/input_files/gwama_input.has${pheno}.topmed74k_topImpukb200k250k_mgbb53k_biovu54k.Aug2023.txt; done < <(echo -e "CH\nDNMT3A\nTET2")

# is file available:  while read pheno; do while read lines; do ll ${lines}; done </medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/input_files/gwama_input.has${pheno}.topmed74k_topImpukb200k250k_mgbb53k_biovu54k.Aug2023.txt; done < <(echo -e "CH\nDNMT3A\nTET2")

## Run Meta-analyses: 
## while read pheno; do qsub -R y -wd /broad/hptmp/mesbah/dataset/ch_gwas/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=40G -l h_rt=10:00:00 -N metaGWAS.has${pheno}.Aug2023.Topmed_UKB_MGB_BioVU /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/01.GWAMA.Meta_Analysis.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/input_files/gwama_input.has${pheno}.topmed74k_topImpukb200k250k_mgbb53k_biovu54k.Aug2023.txt /broad/hptmp/mesbah/dataset/ch_gwas/metaGWAS.has${pheno}.topmed74k_topImpukb200k250k_mgbb53k_biovu54k.Aug2023; done < <(echo -e "CH\nDNMT3A\nTET2")
###
##########################################################################################################


########################################### for CKD MR #####################################
######## August 17, 2023 #########
### Meta-analyse ukb200k and ukb250k GWAS for CKD MR [caitlyn.vlasschaert et al.]
## mkdir -p /broad/hptmp/mesbah/dataset/ch_gwas/ukb_4_CatV
# while read CH_type; do echo -e "/medpop/esp/mesbah/GWAS_CHIP/inputUKB/2023_ukbb/rap_gwas/has${CH_type}.chr1_22.ukb200k.regenie.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/inputUKB/2023_ukbb/rap_gwas/has${CH_type}.chr1_22.ukb250k.regenie.tsv.gz" > /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/input_files/gwama_input.has${CH_type}.topImpUkb200k250k.txt; done < <(echo -e "CH\nDNMT3A\nTET2")

## qsub: while read CH_type; do qsub -R y -wd /broad/hptmp/mesbah/dataset/ch_gwas/ukb_4_CatV -pe smp 1 -binding linear:1 -l h_vmem=40G -l h_rt=10:00:00 -N metaGWAS.${CH_type}.ukb200k250k_TOPMed_imputed /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/01.GWAMA.Meta_Analysis.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/input_files/gwama_input.has${CH_type}.topImpUkb200k250k.txt /broad/hptmp/mesbah/dataset/ch_gwas/ukb_4_CatV/metagwas.has${CH_type}.ukb200k250k_topImp; done < <(echo -e "CH\nDNMT3A\nTET2")

## When done: ####### Prep summary stats 
#### while read pheno; do zcat /broad/hptmp/mesbah/dataset/ch_gwas/ukb_4_CatV/metagwas.has${pheno}.ukb200k250k_topImp.out.gz | awk '(NR==1){print "varID\tEffectAllele\tOtherAllele\tEffectAlleleFrequency\tBETA\tSE\tPval\tN\tEffect_Direction\tQ_Pval\tI_sqr"};(NR>1 && $1~/^chr[0-9]+/ && $4>=0.005 && $4<=0.995 && $15==2){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$10"\t"$16"\t"$17"\t"$13"\t"$14}' | gzip -c >> /broad/hptmp/mesbah/dataset/ch_gwas/ukb_4_CatV/has${pheno}.multiAncestry_metagwas.ukb200k250k.topmedImp.maf005.17Aug2023.tsv.gz; done < <(echo -e "CH\nDNMT3A\nTET2") &
######

############################################################################################

############################### topmed74k + topmed imputed ukb200k + ukb250k + AoU 80k + mgbb53k + bioVu54k
### while read pheno; do echo -e "SNPID\tREF\tALT\tP\tAAF\tBETA\tSE\tN" | gzip -c > /medpop/esp/mesbah/GWAS_CHIP/AoU/chr1_22.AoU_multiAnc_2023.has${pheno}.regenie.tsv.gz; for chr in {1..22}; do bzcat /medpop/esp2/mesbah/projects/Meta_GWAS/2023_gwas/aou_mmu/chr${chr}.AoU_multiAnc_2023.has${pheno}.regenie.tsv.bz2 | awk 'NR>1{print $1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$11}' | sed 's:\_:\::g' | gzip -c >>  /medpop/esp/mesbah/GWAS_CHIP/AoU/chr1_22.AoU_multiAnc_2023.has${pheno}.regenie.tsv.gz; done; done < <(echo -e "CH\nDNMT3A\nTET2") &

## Overall CH
# echo -e "/medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2021/topmed2021_hasCH.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/inputUKB/2023_ukbb/rap_gwas/hasCH.chr1_22.ukb200k.regenie.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/inputUKB/2023_ukbb/rap_gwas/hasCH.chr1_22.ukb250k.regenie.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/AoU/chr1_22.AoU_multiAnc_2023.hasCH.regenie.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/MGBB53k/mgbb53k_gwama.chr1_22.hasCH.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_hasCH_results_merged_subset.tsv.gz" > /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/input_files/gwama_input.hasCH.topmed74k_topImpukb200k250k_aou98k_mgbb53k_biovu54k.txt

## DNMT3A, TET2 
# while read pheno; do echo -e "/medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2021/topmed2021_has${pheno}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/inputUKB/2023_ukbb/rap_gwas/has${pheno}.chr1_22.ukb200k.regenie.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/inputUKB/2023_ukbb/rap_gwas/has${pheno}.chr1_22.ukb250k.regenie.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/AoU/chr1_22.AoU_multiAnc_2023.has${pheno}.regenie.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/MGBB53k/mgbb53k_gwama.chr1_22.has${pheno}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_has${pheno}_results_merged_subset.tsv.gz" > /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/input_files/gwama_input.has${pheno}.topmed74k_topImpukb200k250k_aou98k_mgbb53k_biovu54k.txt; done < <(echo -e "DNMT3A\nTET2")

# while read lines; do ll ${lines}; done </medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/input_files/gwama_input.hasDNMT3A.topmed74k_topImpukb200k250k_aou98k_mgbb53k_biovu54k.txt

# mkdir -p /broad/hptmp/mesbah/ch_gwas/tmpdir

# while read pheno; do qsub -R y -wd /broad/hptmp/mesbah/dataset/ch_gwas/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=40G -l h_rt=10:00:00 -N metaGWAS.${pheno}.T74topImpUK200k250kAoU98kMGB53BVU54 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/01.GWAMA.Meta_Analysis.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/input_files/gwama_input.has${pheno}.topmed74k_topImpukb200k250k_aou98k_mgbb53k_biovu54k.txt /broad/hptmp/mesbah/dataset/ch_gwas/metaGWAS.has${pheno}.topmed74k_topImpukb200k250k_aou98k_mgbb53k_biovu54k; done < <(echo -e "CH\nDNMT3A\nTET2")


###############################

############## Multi-ancestry

##### topmed74k, ukb200k, ukb250k, aou98k, mgbb53k, biovu54k

# while read traits; do echo -e "/medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2021/topmed2021_has${traits}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/inputUKB/allUKB/lifted_hg38.chr1_22.has${traits}.11Aug2021_ukb200k_allsamples.INFO.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/inputUKB/ukb450k/250k/hg38/lifted_hg38.chr1_22.has${traits}.ukb250k.INFO.tsv.gz\n/medpop/esp2/mesbah/datasets/allofus/gwas_ch/${traits}_ch.chr1_22.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/MGBB53k/mgbb53k_gwama.chr1_22.has${traits}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_has${traits}_results_merged_subset.tsv.gz" > /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/input_files/gwama_input.has${traits}.topmed74k_ukb200k250k_aou98k_mgbb53k_biovu54k.txt; done < <(echo -e "CHIP\nDNMT3A\nTET2")

# mkdir -p /broad/hptmp/mesbah/ch_gwas/tmpdir

# while read traits; do qsub -R y -wd /broad/hptmp/mesbah/ch_gwas/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=40G -l h_rt=10:00:00 -N metaGWAS.${traits}.T74UK200k250kAoU98kMGB53BVU54 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/01.GWAMA.Meta_Analysis.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/input_files/gwama_input.has${traits}.topmed74k_ukb200k250k_aou98k_mgbb53k_biovu54k.txt /broad/hptmp/mesbah/ch_gwas/metaGWAS.${traits}.topmed74k_ukb200k250k_aou98k_mgbb53k_biovu54k; done < <(echo -e "CHIP\nDNMT3A\nTET2")


###############################################

## EUR UKB 200k, 250k GWAMA for finemapping
# while read traits; do qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/eur_only_summary -pe smp 1 -binding linear:1 -l h_vmem=30G -l h_rt=10:00:00 -N ukb200k_250k.eur_metaGWAS.${traits} /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/01.GWAMA.Meta_Analysis.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/input_files/ukb200k_ukb250k.meta_${traits}.EUR.txt /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/eur_only_summary/ukb200k_ukb250k.meta_${traits}.EUR.gwama /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/eur_meta_summary.ukb200k_ukb250k.meta_${traits}.tsv; done < <(echo -e "CHIP\nDNMT3A\nTET2")
###
## EUR GWAS
# echo -e "/medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/eur_topmed/topmed2019_EURonly_samples.hasCHIP.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/lifted_hg38.chr1_22.hasCHIP.21Aug2021.UKB200k_onlyEURsamples.INFO.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/inputUKB/ukb450k/250k/hg38/lifted_hg38.chr1_22.hasCHIP.ukb250k_WB.INFO.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/MGBB53k/chr1_22.MGBB53k_EUR.hasCHIP.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_hasCHIP_results_merged_subset.tsv.gz" > /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/input_files/gwama_input.hasCHIP.EUR_only.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.txt
# while read traits; do echo -e "/medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/lifted_hg38.chr1_22.has${traits}.21Aug2021.UKB200k_onlyEURsamples.INFO.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/inputUKB/ukb450k/250k/hg38/lifted_hg38.chr1_22.has${traits}.ukb250k_WB.INFO.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/MGBB53k/chr1_22.MGBB53k_EUR.has${traits}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_has${traits}_results_merged_subset.tsv.gz" > /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/input_files/gwama_input.has${traits}.EUR_only.UKB200kUKB250kMGBB53kBioVU54k.txt; done < <(echo -e "DNMT3A\nTET2")

	# Submit Qsub
# while read traits; do qsub -R y -wd /broad/hptmp/mesbah/gwas/ukb450k/metagwas/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=40G -l h_rt=10:00:00 -N eur_metaGWAS.${traits}.ukb_mgb_biovu /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/01.GWAMA.Meta_Analysis.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/input_files/gwama_input.has${traits}.EUR_only.UKB200kUKB250kMGBB53kBioVU54k.txt /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/eur_metaGWAS.${traits}.ukb_mgb_biovu /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/eur_meta_summary.${traits}.ukb_mgb_biovu.tsv; done < <(echo -e "DNMT3A\nTET2")
# while read traits; do qsub -R y -wd /broad/hptmp/mesbah/gwas/ukb450k/metagwas/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=40G -l h_rt=10:00:00 -N eur_metaGWAS.${traits}.topmed_ukb_mgb_biovu /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/01.GWAMA.Meta_Analysis.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/input_files/gwama_input.hasCHIP.EUR_only.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.txt /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/eur_metaGWAS.${traits}.topmed_ukb_mgb_biovu /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/eur_meta_summary.${traits}.topmed_ukb_mgb_biovu.tsv; done < <(echo -e "CHIP")

################### Multi-Ancestry
### Aug 8, 2022
# TOPMed 2019 ukb200k ukb250k mgbb53k biovu54k
# while read traits; do echo -e "/medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/topmed2019.all_samples.has${traits}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/inputUKB/allUKB/lifted_hg38.chr1_22.has${traits}.11Aug2021_ukb200k_allsamples.INFO.tsv.gz\n/broad/hptmp/mesbah/gwas/ukb450k/250k.step2/ukb250k_gwama_input/lifted_hg38.chr1_22.has${traits}.ukb250k.INFO.tsv.gz\n/broad/hptmp/mesbah/gwas/mgbb53k/mgbb53k_gwama.chr1_22.has${traits}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_has${traits}_results_merged_subset.tsv.gz" > /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/01.PrepSummary/gwama_input.has${traits}.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.txt; done < <(echo -e "CHIP\nDNMT3A\nTET2\nASXL1")
# while read traits; do qsub -R y -wd /broad/hptmp/mesbah/gwas/ukb450k/metagwas/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=40G -l h_rt=10:00:00 -N metaGWAS.${traits}.T64UK200_250kMGB53BVU54 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/01.run.GWAMA.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/01.PrepSummary/gwama_input.has${traits}.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.txt /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/metaGWAS.${traits}.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/meta_summary.${traits}.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.tsv; done < <(echo -e "CHIP\nDNMT3A\nTET2\nASXL1")

# TOPMed UKB450k mgbb53k biovu 54k 
#while read traits; do echo -e "/medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/topmed2019.all_samples.has${traits}.tsv.gz\n/broad/hptmp/mesbah/gwas/ukb450k/450k.step2/lifted_hg38.chr1_22.has${traits}.ukb450k.INFO.tsv.gz\n/broad/hptmp/mesbah/gwas/mgbb53k/mgbb53k_gwama.chr1_22.has${traits}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_has${traits}_results_merged_subset.tsv.gz" > /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/01.PrepSummary/gwama_input.has${traits}.TOPMed2019UKB450kMGBB53kBioVU54k.txt; done < <(echo -e "CHIP\nDNMT3A\nTET2\nASXL1")
# while read traits; do qsub -R y -wd /broad/hptmp/mesbah/gwas/ukb450k/metagwas/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=40G -l h_rt=10:00:00 -N metaGWAS.${traits}.T64UK450kMGB53BVU54 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/01.run.GWAMA.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/01.PrepSummary/gwama_input.has${traits}.TOPMed2019UKB450kMGBB53kBioVU54k.txt /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/metaGWAS.${traits}.TOPMed2019UKB450kMGBB53kBioVU54k; done < <(echo -e "CHIP\nDNMT3A\nTET2\nASXL1")

## ##################
## 26 July 2022
# while read traits; do echo -e "/medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/topmed2019.all_samples.has${traits}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/inputUKB/allUKB/lifted_hg38.chr1_22.has${traits}.11Aug2021_ukb200k_allsamples.INFO.tsv.gz\n/broad/hptmp/mesbah/gwas/mgbb53k/mgbb53k_gwama.chr1_22.has${traits}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_has${traits}_results_merged_subset.tsv.gz" > /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/01.PrepSummary/gwama_input.has${traits}.TOPMed2019UKB200kMGBB53kBioVU54k.txt; done < <(echo -e "CHIP\nDNMT3A\nTET2\nASXL1")
## Jobs killed with h_vmem=20G
# while read traits; do qsub -R y -wd /broad/hptmp/mesbah/gwas/mgbb53k/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=40G -l h_rt=10:00:00 -N metaGWAS.${traits}.TOPMed2019UKB200kMGBB53kBioVU54k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/01.run.GWAMA.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/01.PrepSummary/gwama_input.has${traits}.TOPMed2019UKB200kMGBB53kBioVU54k.txt /broad/hptmp/mesbah/gwas/mgbb53k/metaGWAS/metaGWAS.${traits}.TOPMed2019UKB200kMGBB53kBioVU54k; done < <(echo -e "CHIP\nDNMT3A\nTET2")
##########################

############# old
## EUR 3: while read genename; do qsub -R y -wd /broad/hptmp/mesbah/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=20G -l h_rt=10:00:00 -N gwama_EUR3Samp.${genename} /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_GWAMA/01.run.GWAMA.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_GWAMA/Input_GWAMA.euronly_UKB_MGBB_BioVU.${genename}_GWAS.txt /broad/hptmp/mesbah/ukb_chip/meta_gwas/eurGWAS_meta4/GWAMA.meta_${genename}.euronly_UKB_MGBB_BioVU; done < <(echo -e "CHIP\nDNMT3A\nTET2")

## qsub -R y -wd /broad/hptmp/mesbah/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=30G -l h_rt=10:00:00 -N gwama_allEUR.chip /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_GWAMA/01.run.GWAMA.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_GWAMA/Input_GWAMA.EUR_only_CHIP_GWAS.txt /broad/hptmp/mesbah/ukb_chip/meta_gwas/eurGWAS_meta4/GWAMA.EUR_only_CHIP_GWAS.TopMedUKbbMGbbBioVU
## All: while read genename; do qsub -R y -wd /broad/hptmp/mesbah/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=20G -l h_rt=10:00:00 -N gwama_allSamples.${genename} /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_GWAMA/01.run.GWAMA.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_GWAMA/Input_GWAMA.allsamples.topmed_ukb_mgb_biovu.${genename}_GWAS.txt /broad/hptmp/mesbah/ukb_chip/meta_gwas/meta4.topukbmgbbiovuSaige/GWAMA.allsamples.meta_${genename}.TopMedUKbbMGbbBioVU; done < <(echo -e "CHIP\nDNMT3A\nTET2")

## Mixed: all + EUR
# while read genename; do qsub -R y -wd /broad/hptmp/mesbah/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=20G -l h_rt=10:00:00 -N gwama_mixedSamples.${genename} /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_GWAMA/01.run.GWAMA.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_GWAMA/Input_GWAMA.allTOPMed_eurUKB_MGBB_BioVU.${genename}_GWAS.txt /broad/hptmp/mesbah/ukb_chip/meta_gwas/eurGWAS_meta4/GWAMA.mixedSamples.meta_${genename}.allTOPmed_EUR_UKB_MGB_BioVU; done < <(echo -e "CHIP\nDNMT3A\nTET2")

## 
gwas_list=${1}

output_prefix=${2}

# meta_sum=${3} # summary tsv file for plotting

################
## Prepare GWAMA input
# mkdir -p /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/{aou,topmed,ukbb,mgbb,biovu}
# BioVU: 
# while read pheno; do zcat /medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_has${pheno}_results_merged_subset.tsv.gz | awk -v OFS='\t' '{print $3,$4,$5,$7,$10,$11,$9,$13}' > /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/biovu/BioVU.saige_has${pheno}_results_merged_subset.tsv; done < <(echo -e "CH\nDNMT3A\nTET2") &
# MGBB
# while read ANC; do while read pheno; do zcat /medpop/esp/mesbah/GWAS_CHIP/MGBB53k/${ANC}/has${pheno}.chr1_22.${ANC}.regenie.tsv.gz > /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/mgbb/has${pheno}.chr1_22.${ANC}.regenie.tsv; done < <(echo -e "CH\nDNMT3A\nTET2\nCHvaf10\nASXL1\nSF\nDDR"); done < <(echo -e "MultiANC\nFemale\nMale\nEUR\nAFR\nAMR") &
# UKBB

# TOPMed
# while read ANC; do while read pheno; do zcat /broad/hptmp/mesbah/dataset/ch_gwas/topmed/gwas/2024/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv.gz > /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.tsv; done < <(echo -e "CH\nDNMT3A\nTET2\nCHvaf10\nASXL1\nSF\nDDR"); done < <(echo -e "MultiANC\nFemale\nMale\nEUR\nAFR\nAMR") &
# AoU


#
###
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
bgzip -f ${output_prefix}.out

###### Meta-Analysis Summary
# min MAF >=0.001 && MAF <= 0.999 && N_Studies>1
### Extart SNP present in >=2 studies, MAF>=0.1%
# zcat  metaGWAS.CHIP.TOPMed2019UKB200kMGBB53kBioVU54k.out.gz | awk '(NR==1){print "RSID\tREF\tALT\tAAF\tBETA\tSE\tP"}(NR>1 && $1~/^chr[0-9]+/ && $15>1 && $4>=0.001 && $4<=0.999){print $1"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6"\t"$10}' | bgzip -c > chr1_22.n2maf001.CHIP.tsv.gz &
# zcat  metaGWAS.DNMT3A.TOPMed2019UKB200kMGBB53kBioVU54k.out.gz | awk '(NR==1){print "RSID\tREF\tALT\tAAF\tBETA\tSE\tP"}(NR>1 && $1~/^chr[0-9]+/ && $15>1 && $4>=0.001 && $4<=0.999){print $1"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6"\t"$10}' | bgzip -c > chr1_22.n2maf001.DNMT3A.tsv.gz &
# zcat  metaGWAS.TET2.TOPMed2019UKB200kMGBB53kBioVU54k.out.gz | awk '(NR==1){print "RSID\tREF\tALT\tAAF\tBETA\tSE\tP"}(NR>1 && $1~/^chr[0-9]+/ && $15>1 && $4>=0.001 && $4<=0.999){print $1"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6"\t"$10}' | bgzip -c > chr1_22.n2maf001.TET2.tsv.gz &
# zcat chr1_22.n2maf001.CHIP.tsv.gz | head -1 > sorted.chr1_22.n2maf001.CHIP.tsv && zcat chr1_22.n2maf001.CHIP.tsv.gz | awk 'NR>1' | sort -k1 -V >> sorted.chr1_22.n2maf001.CHIP.tsv && bgzip sorted.chr1_22.n2maf001.CHIP.tsv &
####### 

##### Uncomment to run following code
# echo -e "CHR\tPOS\tREF\tALT\tvarID\tAAF\tBETA\tSE\tP" > ${meta_sum}

# zcat ${output_prefix}.out.gz | awk '(NR>1 && $1~/^chr[0-9]+/ && $15>1 && $4>=0.001 && $4<=0.999){print $1"\t"$4"\t"$5"\t"$6"\t"$10}' | sort -k1 -V | sed 's:\::\t:g' | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$1":"$2":"$3":"$4"\t"$5"\t"$6"\t"$7"\t"$8}' >> ${meta_sum}

# bgzip -f ${meta_sum} 

# tabix -f -s 1 -b 2 -e 2 ${meta_sum}.gz

###########


######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

