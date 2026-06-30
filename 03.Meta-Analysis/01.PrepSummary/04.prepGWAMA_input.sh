#!/bin/bash

## qsub -R y -wd /broad/hptmp/mesbah/tmpdir -N prep_GWAMA -pe smp 1 -binding linear:1 -l h_rt=20:00:00 -l h_vmem=30G /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/01.PrepSummary/04.prepGWAMA_input.sh

source /broad/software/scripts/useuse

use Tabix

######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
#########
	
	## TOPMED
# echo -e "CHR\tPOS\tSNPID\tREF\tALT\tAC_ALT\tAAF\tN\tBETA\tSE\tTstat\tP\tp.value.NA\tIs.SPA.converge\tvarT\tvarTstar" > /medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/eur_topmed/topmed2019_EURonly_samples.hasCHIP.tsv
# zcat /medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/eur_topmed/topmed2019_EUR.hasCHIP.tsv.gz | awk 'NR>1' >> /medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/eur_topmed/topmed2019_EURonly_samples.hasCHIP.tsv
# bgzip -f /medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/eur_topmed/topmed2019_EURonly_samples.hasCHIP.tsv
# tabix -f -s 1 -b 2 -e 2 /medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/eur_topmed/topmed2019_EURonly_samples.hasCHIP.tsv.gz
##
# while read phenotype
# do
# 	echo -e "CHR\tPOS\tSNPID\tREF\tALT\tAC_ALT\tAAF\tN\tBETA\tSE\tTstat\tP\tp.value.NA\tIs.SPA.converge\tvarT\tvarTstar" > /medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/topmed2019.all_samples.has${phenotype}.tsv
# 	zcat /medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/topmed2019.has${phenotype}.tsv.gz | awk 'NR>1' >> /medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/topmed2019.all_samples.has${phenotype}.tsv
# 	bgzip -f /medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/topmed2019.all_samples.has${phenotype}.tsv
# 	tabix -f -s 1 -b 2 -e 2 /medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/topmed2019.all_samples.has${phenotype}.tsv.gz 
# done < <(echo -e "CHIP\nDNMT3A\nTET2" )

	## UKBB
while read phenotype 
do
	echo -e "SNPID\tCHR\tPOS\tREF\tALT\tAAF\tBETA\tSE\tP\tN\told_hg19_id" > /medpop/esp/mesbah/GWAS_CHIP/inputUKB/allUKB/lifted_hg38.chr1_22.has${phenotype}.11Aug2021_ukb200k_allsamples.INFO.tsv

	zcat /medpop/esp/mesbah/GWAS_CHIP/inputUKB/allUKB/lifted_hg38.chr1_22.has${phenotype}.11Aug2021_ukb200k.INFO.tsv.gz | awk 'NR>1' >> /medpop/esp/mesbah/GWAS_CHIP/inputUKB/allUKB/lifted_hg38.chr1_22.has${phenotype}.11Aug2021_ukb200k_allsamples.INFO.tsv

	bgzip -f /medpop/esp/mesbah/GWAS_CHIP/inputUKB/allUKB/lifted_hg38.chr1_22.has${phenotype}.11Aug2021_ukb200k_allsamples.INFO.tsv
	
	tabix -f -s 2 -b 3 -e 3 /medpop/esp/mesbah/GWAS_CHIP/inputUKB/allUKB/lifted_hg38.chr1_22.has${phenotype}.11Aug2021_ukb200k_allsamples.INFO.tsv.gz
	
	# EUR
	echo -e "SNPID\tCHR\tPOS\tREF\tALT\tAAF\tBETA\tSE\tP\tN\told_hg19_id" > /medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/lifted_hg38.chr1_22.has${phenotype}.21Aug2021.UKB200k_onlyEURsamples.INFO.tsv

	zcat /medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/lifted_hg38.chr1_22.has${phenotype}.21Aug2021_ukbEUR.INFO.tsv.gz | awk 'NR>1' >> /medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/lifted_hg38.chr1_22.has${phenotype}.21Aug2021.UKB200k_onlyEURsamples.INFO.tsv
	
	bgzip -f /medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/lifted_hg38.chr1_22.has${phenotype}.21Aug2021.UKB200k_onlyEURsamples.INFO.tsv
	
	tabix -f -s 2 -b 3 -e 3 /medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/lifted_hg38.chr1_22.has${phenotype}.21Aug2021.UKB200k_onlyEURsamples.INFO.tsv.gz

done < <(echo -e "CHIP\nDNMT3A\nTET2")


	## MGBB
while read phenotype 
do 
	echo -e "CHR\tPOS\tSNPID\tREF\tALT\tTrait\tCohort\tModel\tEffect\tLCI_Effect\tUCI_Effect\tP\tAAF\tNum_Cases\tCases_Ref\tCases_Het\tCases_Alt\tNum_Controls\tControls_Ref\tControls_Het\tControls_Alt\tBETA\tSE\tMAC\tN" > /medpop/esp/mesbah/GWAS_CHIP/MGB_GSA/allMGBB/MGB_all_samples.chr1_22.has${phenotype}.tsv
	zcat /medpop/esp/mesbah/GWAS_CHIP/MGB_GSA/allMGBB/MGB_all.chr1_22.has${phenotype}.tsv.gz | awk 'NR>1' >> /medpop/esp/mesbah/GWAS_CHIP/MGB_GSA/allMGBB/MGB_all_samples.chr1_22.has${phenotype}.tsv
	bgzip -f /medpop/esp/mesbah/GWAS_CHIP/MGB_GSA/allMGBB/MGB_all_samples.chr1_22.has${phenotype}.tsv
	tabix -f -s 1 -b 2 -e 2 /medpop/esp/mesbah/GWAS_CHIP/MGB_GSA/allMGBB/MGB_all_samples.chr1_22.has${phenotype}.tsv.gz
	# EUR samples
	echo -e "CHR\tPOS\tSNPID\tREF\tALT\tTrait\tCohort\tModel\tEffect\tLCI_Effect\tUCI_Effect\tP\tAAF\tNum_Cases\tCases_Ref\tCases_Het\tCases_Alt\tNum_Controls\tControls_Ref\tControls_Het\tControls_Alt\tBETA\tSE\tMAC\tN" > /medpop/esp/mesbah/GWAS_CHIP/MGB_GSA/eurMGBB/MGB_onlyEURsamples.chr1_22.has${phenotype}.tsv
	zcat /medpop/esp/mesbah/GWAS_CHIP/MGB_GSA/eurMGBB/MGB_EUR.chr1_22.has${phenotype}.tsv.gz | awk 'NR>1' >> /medpop/esp/mesbah/GWAS_CHIP/MGB_GSA/eurMGBB/MGB_onlyEURsamples.chr1_22.has${phenotype}.tsv
	bgzip -f /medpop/esp/mesbah/GWAS_CHIP/MGB_GSA/eurMGBB/MGB_onlyEURsamples.chr1_22.has${phenotype}.tsv
	tabix -f -s 1 -b 2 -e 2 /medpop/esp/mesbah/GWAS_CHIP/MGB_GSA/eurMGBB/MGB_onlyEURsamples.chr1_22.has${phenotype}.tsv.gz
done < <(echo -e "CHIP\nDNMT3A\nTET2")

	## BioVU
while read phenotype 
do 
	echo -e "CHR\tPOS\tSNPID\tREF\tALT\tAC_ALT\tAAF\timputationInfo\tN\tBETA\tSE\tTstat\tP\tp.value.NA\tIs.SPA.converge\tvarvarTstar" > /medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_has${phenotype}_results_merged_subset.tsv
	
	zcat /medpop/esp/mesbah/GWAS_CHIP/BioVU/biovu.saige_has${phenotype}_results_merged_subset.tsv.gz | awk 'NR>1' >> /medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_has${phenotype}_results_merged_subset.tsv
	
	bgzip -f /medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_has${phenotype}_results_merged_subset.tsv
	
	tabix -f -s 1 -b 2 -e 2 /medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_has${phenotype}_results_merged_subset.tsv.gz
done < <(echo -e "CHIP\nDNMT3A\nTET2")

##### Prepare INPUT files

echo -e "/medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/eur_topmed/topmed2019_EURonly_samples.hasCHIP.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/lifted_hg38.chr1_22.hasCHIP.21Aug2021.UKB200k_onlyEURsamples.INFO.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/MGB_GSA/eurMGBB/MGB_onlyEURsamples.chr1_22.hasCHIP.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_hasCHIP_results_merged_subset.tsv.gz" > /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_GWAMA/Input_GWAMA.EUR_only_CHIP_GWAS.txt

while read phenotype
do
	echo -e "/medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/topmed2019.all_samples.has${phenotype}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/inputUKB/allUKB/lifted_hg38.chr1_22.has${phenotype}.11Aug2021_ukb200k_allsamples.INFO.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/MGB_GSA/allMGBB/MGB_all_samples.chr1_22.has${phenotype}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_has${phenotype}_results_merged_subset.tsv.gz" > /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_GWAMA/Input_GWAMA.allsamples.topmed_ukb_mgb_biovu.${phenotype}_GWAS.txt
	
	echo -e "/medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/topmed2019.all_samples.has${phenotype}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/lifted_hg38.chr1_22.has${phenotype}.21Aug2021.UKB200k_onlyEURsamples.INFO.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/MGB_GSA/eurMGBB/MGB_onlyEURsamples.chr1_22.has${phenotype}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_has${phenotype}_results_merged_subset.tsv.gz" > /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_GWAMA/Input_GWAMA.allTOPMed_eurUKB_MGBB_BioVU.${phenotype}_GWAS.txt

done < <(echo -e "CHIP\nDNMT3A\nTET2")

# compress 

#####

######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########


