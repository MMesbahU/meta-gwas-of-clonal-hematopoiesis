#!/bin/bash
###
## Oct 2023: AoU250k UKB200k UKB250k TOPMed 74k MGB53k BioVU54k
## Metal path: #/home/jupyter/workspaces/detectionofclonalhematopoiesisofindeterminatepotentialchip/metal/random-metal-0.1.0/executables/metal
#Random metal: /home/jupyter/workspaces/detectionofclonalhematopoiesisofindeterminatepotentialchip/metal/METAL-2020-05-05/metal/bin/metal
## Run: 
# while read phenotype; do bash /home/jupyter/myscripts/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/02.METAL.Meta_Analysis.sh /home/jupyter/myscripts/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/input_files/metal_directive.AoU250k_UKB200k250k_TOPMed74k_MGB53k_BioVU54k.multiAnc.${phenotype}.txt /home/jupyter/ch_gwas/aou250k/aou_gwas_summary/chr1_22.AoU_multiAnc_2023.has${phenotype}.August2023.regenie.tsv.gz /home/jupyter/ch_gwas/aou250k/sum4_cohort/has${phenotype}.chr1_22.ukb200k.regenie.tsv.gz /home/jupyter/ch_gwas/aou250k/sum4_cohort/has${phenotype}.chr1_22.ukb250k.regenie.tsv.gz /home/jupyter/ch_gwas/aou250k/sum4_cohort/topmed2021_has${phenotype}.tsv.gz /home/jupyter/ch_gwas/aou250k/sum4_cohort/mgbb53k_gwama.chr1_22.has${phenotype}.tsv.gz /home/jupyter/ch_gwas/aou250k/sum4_cohort/BioVU.saige_has${phenotype}_results_merged_subset.tsv.gz /home/jupyter/ch_gwas/aou250k/gwama_meta/metal_summary/metalGWAS.multiAnc.${phenotype}.AoU250k_UKB200k250k_TOPMed74k_MGB53k_BioVU54k.Oct2023 /home/jupyter/ch_gwas/aou250k/gwama_meta/metal_summary /home/jupyter/myscripts/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/01.plot.metaGWAS.R /home/jupyter/ch_gwas/aou250k/gwama_meta/metal_summary/summuray_for_plot.AoU250k_UKB200k250k_TOPMed74k_MGB53k_BioVU54k.multiAnc.${phenotype}.tsv ${phenotype} /home/jupyter/ch_gwas/aou250k/gwama_meta/metal_summary/unsorted_metal.AoU250k_UKB200k250k_TOPMed74k_MGB53k_BioVU54k.multiAnc.${phenotype}.vcf /home/jupyter/ch_gwas/aou250k/gwama_meta/metal_summary/sorted_metal.AoU250k_UKB200k250k_TOPMed74k_MGB53k_BioVU54k.multiAnc.${phenotype}.vcf /home/jupyter/workspaces/detectionofclonalhematopoiesisofindeterminatepotentialchip/metal/random-metal-0.1.0/executables/metal; done < <(echo -e "CH\nDNMT3A\nTET2")
###############
#### Multi-ancestry 
## 
## TOPMed2021 UKB200k UKB250 AoU98k MGBB53k BioVU54k
# updated 
## topmed imputed UKB RAP 
## AoU 98k mmu gwas
## Topmed 74k GWAS
###
# while read phenotype; do qsub -R y -wd /broad/hptmp/mesbah/dataset/ch_gwas/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=40G -l h_rt=20:00:00 -N meta_all${phenotype} /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/02.METAL.Meta_Analysis.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/input_files/metal_directive.topmed74k_ukb200k250k_aou98k_mgbb53k_biovu54k.multiAnc.${phenotype}.txt /medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2021/topmed2021_has${phenotype}.tsv.gz /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2023_ukbb/rap_gwas/has${phenotype}.chr1_22.ukb200k.regenie.tsv.gz /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2023_ukbb/rap_gwas/has${phenotype}.chr1_22.ukb250k.regenie.tsv.gz /medpop/esp/mesbah/GWAS_CHIP/AoU/chr1_22.AoU_multiAnc_2023.has${phenotype}.regenie.tsv.gz /medpop/esp/mesbah/GWAS_CHIP/MGBB53k/mgbb53k_gwama.chr1_22.has${phenotype}.tsv.gz /medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_has${phenotype}_results_merged_subset.tsv.gz /broad/hptmp/mesbah/dataset/ch_gwas/metalGWAS.multiAnc.${phenotype}.topmed74k_ukb200k250k_aou98k_mgbb53k_biovu54k /broad/hptmp/mesbah/dataset/ch_gwas /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/01.plot.metaGWAS.R /broad/hptmp/mesbah/dataset/ch_gwas/summuray_for_plot.topmed74k_ukb450k_aou98k_mgbb53k_biovu54k.multiAnc.${phenotype}.tsv ${phenotype} /broad/hptmp/mesbah/dataset/ch_gwas/unsorted_metal.topmed74k_ukb200k250k_aou98k_mgbb53k_biovu54k.multiAnc.${phenotype}.vcf /broad/hptmp/mesbah/dataset/ch_gwas/sorted_metal.topmed74k_ukb200k250k_aou98k_mgbb53k_biovu54k.multiAnc.${phenotype}.vcf; done < <(echo -e "CH\nDNMT3A\nTET2")

##
# Metal + Plot: 
# while read phenotype; do qsub -R y -wd /broad/hptmp/mesbah/ch_gwas/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=40G -l h_rt=20:00:00 -N meta_all${phenotype} /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/02.METAL.Meta_Analysis.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/input_files/metal_directive.topmed74k_ukb200k250k_aou98k_mgbb53k_biovu54k.${phenotype}.txt /medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2021/topmed2021_has${phenotype}.tsv.gz /medpop/esp/mesbah/GWAS_CHIP/inputUKB/allUKB/lifted_hg38.chr1_22.has${phenotype}.11Aug2021_ukb200k_allsamples.INFO.tsv.gz /medpop/esp/mesbah/GWAS_CHIP/inputUKB/ukb450k/250k/hg38/lifted_hg38.chr1_22.has${phenotype}.ukb250k.INFO.tsv.gz /medpop/esp2/mesbah/datasets/allofus/gwas_ch/${phenotype}_ch.chr1_22.tsv.gz /medpop/esp/mesbah/GWAS_CHIP/MGBB53k/mgbb53k_gwama.chr1_22.has${phenotype}.tsv.gz /medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_has${phenotype}_results_merged_subset.tsv.gz /broad/hptmp/mesbah/ch_gwas/metaGWAS.${phenotype}.topmed74k_ukb200k250k_aou98k_mgbb53k_biovu54k /broad/hptmp/mesbah/ch_gwas /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/01.plot.metaGWAS.R /broad/hptmp/mesbah/ch_gwas/summuray_for_plot.topmed74k_ukb450k_aou98k_mgbb53k_biovu54k.${phenotype}.tsv ${phenotype} /broad/hptmp/mesbah/ch_gwas/unsorted_metal.topmed74k_ukb200k250k_aou98k_mgbb53k_biovu54k.${phenotype}.vcf /broad/hptmp/mesbah/ch_gwas/sorted_metal.topmed74k_ukb200k250k_aou98k_mgbb53k_biovu54k.${phenotype}.vcf; done < <(echo -e "CHIP\nDNMT3A\nTET2")
####
### EUR only
## while read phenotype; do qsub -R y -wd /broad/hptmp/mesbah/gwas/ukb450k/metagwas/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=40G -l h_rt=20:00:00 -N meta_eur.${phenotype} /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/02.METAL_all.Meta_Analysis.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/metal_directive_eur.ukb450k_mgbb53k_biovu54k.${phenotype}.txt /medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/eur_topmed/topmed2019_EURonly_samples.has${phenotype}.tsv.gz /medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/lifted_hg38.chr1_22.has${phenotype}.21Aug2021.UKB200k_onlyEURsamples.INFO.tsv.gz /broad/hptmp/mesbah/gwas/ukb450k/250k.step2/ukb250k_WB/lifted_hg38.chr1_22.has${phenotype}.ukb250k_WB.INFO.tsv.gz /medpop/esp/mesbah/GWAS_CHIP/MGBB53k/chr1_22.MGBB53k_EUR.has${phenotype}.tsv.gz /medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_has${phenotype}_results_merged_subset.tsv.gz /broad/hptmp/mesbah/gwas/ukb450k/metagwas/metal_GWAS_eur.${phenotype}.topmed64k_ukb450k_mgbb53k_biovu54k /broad/hptmp/mesbah/gwas/ukb450k/metagwas /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/01.plot.metaGWAS.R /broad/hptmp/mesbah/gwas/ukb450k/metagwas/summuray_for_plot_eur.topmed_ukb450k_mgbb53k_biovu54k.${phenotype}.tsv ${phenotype} /broad/hptmp/mesbah/gwas/ukb450k/metagwas/unsorted_metal_eur.topmed_ukb450k_mgbb53k_biovu54k.${phenotype}.vcf /broad/hptmp/mesbah/gwas/ukb450k/metagwas/sorted_metal_eur.topmed_ukb450k_mgbb53k_biovu54k.${phenotype}.vcf; done < <(echo -e "CHIP\nDNMT3A\nTET2")

## TOPMed 2019 UKB200k UKB250 MGBB53k BioVU54k
# Metal + Plot: while read phenotype; do qsub -R y -wd /broad/hptmp/mesbah/gwas/ukb450k/metagwas/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=40G -l h_rt=20:00:00 -N meta_all${phenotype} /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/02.METAL_all.Meta_Analysis.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/metal_directive.topmed64k_ukb450k_mgbb53k_biovu54k.${phenotype}.txt /medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/topmed2019.all_samples.has${phenotype}.tsv.gz /medpop/esp/mesbah/GWAS_CHIP/inputUKB/allUKB/lifted_hg38.chr1_22.has${phenotype}.11Aug2021_ukb200k_allsamples.INFO.tsv.gz /broad/hptmp/mesbah/gwas/ukb450k/250k.step2/ukb250k_gwama_input/lifted_hg38.chr1_22.has${phenotype}.ukb250k.INFO.tsv.gz /broad/hptmp/mesbah/gwas/mgbb53k/mgbb53k_gwama.chr1_22.has${phenotype}.tsv.gz /medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_has${phenotype}_results_merged_subset.tsv.gz /broad/hptmp/mesbah/gwas/ukb450k/metagwas/metaGWAS.${phenotype}.topmed64k_ukb450k_mgbb53k_biovu54k /broad/hptmp/mesbah/gwas/ukb450k/metagwas /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/01.plot.metaGWAS.R /broad/hptmp/mesbah/gwas/ukb450k/metagwas/summuray_for_plot.topmed64k_ukb450k_mgbb53k_biovu54k.${phenotype}.tsv ${phenotype} /broad/hptmp/mesbah/gwas/ukb450k/metagwas/unsorted_metal.topmed64k_ukb450k_mgbb53k_biovu54k.${phenotype}.vcf /broad/hptmp/mesbah/gwas/ukb450k/metagwas/sorted_metal.topmed64k_ukb450k_mgbb53k_biovu54k.${phenotype}.vcf; done < <(echo -e "CHIP\nDNMT3A\nTET2")
#############################

### All Samples: 
# Metal + Plot: while read phenotype; do qsub -R y -wd /broad/hptmp/mesbah/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=40G -l h_rt=20:00:00 -N meta_all${phenotype} /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/02.METAL_all.Meta_Analysis.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_Metal/metal_directive.allSamples.${phenotype}.txt /medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/topmed2019.has${phenotype}.tsv.gz /medpop/esp/mesbah/GWAS_CHIP/inputUKB/allUKB/lifted_hg38.chr1_22.has${phenotype}.11Aug2021_ukb200k.INFO.tsv.gz /medpop/esp/mesbah/GWAS_CHIP/MGB_GSA/allMGBB/MGB_all.chr1_22.has${phenotype}.tsv.gz /medpop/esp/mesbah/GWAS_CHIP/BioVU/biovu.saige_has${phenotype}_results_merged_subset.tsv.gz /broad/hptmp/mesbah/ukb_chip/meta_gwas/meta4.topukbmgbbiovuSaige/meta4.${phenotype}.allSamples_topmed_ukb_mgbb_biovu /broad/hptmp/mesbah/ukb_chip/meta_gwas/meta4.topukbmgbbiovuSaige /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/02.plotGWAS/01.plot.metaGWAS.R /broad/hptmp/mesbah/ukb_chip/meta_gwas/meta4.topukbmgbbiovuSaige/summuray_for_plot.allSamples${phenotype}.tsv ${phenotype}; done < <(echo -e "CHIP\nDNMT3A\nTET2")

## All EUR + TOPMED all: while read phenotype; do qsub -R y -wd /broad/hptmp/mesbah/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=40G -l h_rt=20:00:00 -N meta_alltopmed_3EUR${phenotype} /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/02.METAL_all.Meta_Analysis.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_Metal/metal_directive.3EURplusallTOPMEDsamples.${phenotype}.txt /medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/topmed2019.has${phenotype}.tsv.gz /medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/lifted_hg38.chr1_22.has${phenotype}.21Aug2021_ukbEUR.INFO.tsv.gz /medpop/esp/mesbah/GWAS_CHIP/MGB_GSA/eurMGBB/MGB_EUR.chr1_22.has${phenotype}.tsv.gz /medpop/esp/mesbah/GWAS_CHIP/BioVU/biovu.saige_has${phenotype}_results_merged_subset.tsv.gz /broad/hptmp/mesbah/ukb_chip/meta_gwas/eurGWAS_meta4/meta4.${phenotype}.alltopmed_3EUR.ukb_mgbb_biovu /broad/hptmp/mesbah/ukb_chip/meta_gwas/eurGWAS_meta4 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/02.plotGWAS/01.plot.metaGWAS.R /broad/hptmp/mesbah/ukb_chip/meta_gwas/eurGWAS_meta4/summuray_for_plot.alltopmed_3EUR.${phenotype}.tsv ${phenotype}; done < <(echo -e "CHIP\nDNMT3A\nTET2")

## ALL EUR for CHIP only: phenotype="CHIP"; qsub -R y -wd /broad/hptmp/mesbah/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=30G -l h_rt=10:00:00 -N meta_allEUR.${phenotype} /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/02.METAL_all.Meta_Analysis.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_Metal/metal_directive.allEURsamples.${phenotype}.txt /medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/eur_topmed/topmed2019_EUR.has${phenotype}.tsv.gz /medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/lifted_hg38.chr1_22.has${phenotype}.21Aug2021_ukbEUR.INFO.tsv.gz /medpop/esp/mesbah/GWAS_CHIP/MGB_GSA/eurMGBB/MGB_EUR.chr1_22.has${phenotype}.tsv.gz /medpop/esp/mesbah/GWAS_CHIP/BioVU/biovu.saige_has${phenotype}_results_merged_subset.tsv.gz /broad/hptmp/mesbah/ukb_chip/meta_gwas/eurGWAS_meta4/meta4.${phenotype}.allEUR.topmed_ukb_mgbb_biovu /broad/hptmp/mesbah/ukb_chip/meta_gwas /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/02.plotGWAS/01.plot.metaGWAS.R /broad/hptmp/mesbah/ukb_chip/meta_gwas/eurGWAS_meta4/summuray_for_plot.allEUR.${phenotype}.tsv ${phenotype}

########################

# source /broad/software/scripts/useuse

# use Tabix

# use VCFtools

# use R-4.0

# 
metal_directive=${1}
AoU=${2}
ukb200k=${3}
ukb250k=${4}
topmed=${5}
MGBB=${6}
BioVU=${7}
output_meta=${8}

## Plot input
plotOutdir=${9}
plot_script=${10}
summary_for_plot=${11}
myTrait=${12}

 ## VCF input
unsorted_Metal_vcf=${13}
sorted_Metal_vcf=${14}
# prep_meta_vcf=${12}
Metal=${15}
######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
#########

####
# Metal=/medpop/esp2/mesbah/tools/METAL/build/bin/metal
# Metal=/medpop/esp2/mesbah/tools/random-metal/executables/metal

################################################## Metal ######################
echo "# Meta-analysis weighted by standard error" > ${metal_directive}
echo "SCHEME	STDERR" >> ${metal_directive}

echo "# on genomic control" >> ${metal_directive}
echo "# GENOMICCONTROL ON" >> ${metal_directive}

echo "# allele frequencies in the meta-analysis." >> ${metal_directive}
echo "AVERAGEFREQ ON" >> ${metal_directive}
echo "MINMAXFREQ ON" >> ${metal_directive}

echo "## Track Total N" >> ${metal_directive}
echo "CUSTOMVARIABLE TotalSampleSize" >> ${metal_directive}

echo "####################################################### All of Us 250k WGS GWAS Summary ######################" >> ${metal_directive}
##### Describe and process the AoU input files
#######################################
echo "MARKER   SNPID" >> ${metal_directive}
echo "WEIGHT   N" >> ${metal_directive}
echo "ALLELE   ALT REF" >> ${metal_directive}
echo "FREQ     AAF" >> ${metal_directive}
echo "EFFECT   BETA" >> ${metal_directive}
echo "STDERR   SE" >> ${metal_directive}
# echo "PVAL     P" >> ${metal_directive}
echo "SEPARATOR   TAB" >> ${metal_directive}
echo "LABEL TotalSampleSize as N" >> ${metal_directive}
echo "PROCESS ${AoU}" >> ${metal_directive}
echo "#####################################################################################################" >> ${metal_directive}


echo "############################################## UKB 200k File ##############################" >> ${metal_directive}
echo "SEPARATOR   TAB" >> ${metal_directive}
echo "MARKER   SNPID" >> ${metal_directive}
echo "ALLELE   ALT REF" >> ${metal_directive}
echo "FREQ     AAF" >> ${metal_directive}
echo "WEIGHT N" >> ${metal_directive}
echo "EFFECT   BETA" >> ${metal_directive}
echo "STDERR   SE" >> ${metal_directive}
# echo "PVAL    PVAL" >> ${metal_directive}
echo "LABEL TotalSampleSize as N" >> ${metal_directive}
echo "PROCESS ${ukb200k}" >> ${metal_directive}
echo "############################################################################################" >> ${metal_directive}

echo "############################################## UKB 250k File ##############################" >> ${metal_directive}
echo "SEPARATOR   TAB" >> ${metal_directive}
echo "MARKER   SNPID" >> ${metal_directive}
echo "ALLELE   ALT REF" >> ${metal_directive}
echo "FREQ     AAF" >> ${metal_directive}
echo "WEIGHT N" >> ${metal_directive}
echo "EFFECT   BETA" >> ${metal_directive}
echo "STDERR   SE" >> ${metal_directive}
# echo "PVAL    PVAL" >> ${metal_directive}
echo "LABEL TotalSampleSize as N" >> ${metal_directive}
echo "PROCESS ${ukb250k}" >> ${metal_directive}
echo "############################################################################################" >> ${metal_directive}

echo "############################################## TOPMed 28k File ##############################" >> ${metal_directive}
echo "SEPARATOR   TAB" >> ${metal_directive}
echo "MARKER   SNPID" >> ${metal_directive}
echo "ALLELE   ALT REF" >> ${metal_directive}
echo "FREQ     AAF" >> ${metal_directive}
echo "WEIGHT   N" >> ${metal_directive}
echo "EFFECT   BETA" >> ${metal_directive}
echo "STDERR   SE" >> ${metal_directive}
# echo "PVAL    p.value" >> ${metal_directive}
echo "LABEL TotalSampleSize as N" >> ${metal_directive}
echo "PROCESS ${topmed}" >> ${metal_directive}
echo "############################################################################################" >> ${metal_directive}


echo "####################################################### MGB GSA 53k GWAS File ######################" >> ${metal_directive}
##### Describe and process the MGB input files
#######################################
echo "MARKER   SNPID" >> ${metal_directive}
echo "WEIGHT   N" >> ${metal_directive}
echo "ALLELE   ALT REF" >> ${metal_directive}
echo "FREQ     AAF" >> ${metal_directive}
echo "EFFECT   BETA" >> ${metal_directive}
echo "STDERR   SE" >> ${metal_directive}
# echo "PVAL     Pval" >> ${metal_directive}
echo "SEPARATOR   TAB" >> ${metal_directive}
echo "LABEL TotalSampleSize as N" >> ${metal_directive}
echo "PROCESS ${MGBB}" >> ${metal_directive}
echo "#####################################################################################################" >> ${metal_directive}

echo "############################### BioVU #####################################" >> ${metal_directive}
echo "### Contact: Bick, Alexander <alexander.bick@vumc.org> and Kresge, Hailey A <hailey.a.kresge@vanderbilt.edu>" >> ${metal_directive}
echo "MARKER   SNPID" >> ${metal_directive}
echo "WEIGHT   N" >> ${metal_directive}
echo "ALLELE   ALT REF" >> ${metal_directive}
echo "FREQ     AAF" >> ${metal_directive}
echo "EFFECT   BETA" >> ${metal_directive}
echo "STDERR   SE" >> ${metal_directive}
# echo "PVAL     p.value" >> ${metal_directive}
echo "SEPARATOR   TAB" >> ${metal_directive}
echo "LABEL TotalSampleSize as N" >> ${metal_directive}
echo "PROCESS ${BioVU}" >> ${metal_directive}
echo "##################################################################################" >> ${metal_directive}
echo "# Execute meta-analysis" >> ${metal_directive}

echo "OUTFILE ${output_meta} .tsv" >> ${metal_directive}

echo "ANALYZE HETEROGENEITY" >> ${metal_directive}

# echo "ANALYZE RANDOM" >> ${metal_directive}

echo "QUIT" >> ${metal_directive}

#############################################################

######################### Run Metal
${Metal} ${metal_directive}
########################

######## Plot Manhattan ###################################
# prep for plot
# cut -f1,4,10 ${output_meta}1.tsv | awk 'NR>1 && $2>0.001 && $2<0.999' | sed 's:\::\t:g' | awk '{print $1":"$2":"$3":"$4"\t"$1"\t"$2"\t"$6}' > ${summary_for_plot}

### Run Rscript
# Rscript ${plot_script} ${summary_for_plot} ${plotOutdir} ${myTrait}
##########################################################


## compress
# gzip -f ${output_meta}1.tsv
# gzip -f ${summary_for_plot}
###


########################## Make Metal output VCF
	# VCF Header
echo -e '##fileformat=VCFv4.2' > ${unsorted_Metal_vcf}
echo -e '##INFO=<ID=P_METAL,Number=A,Type=Float,Description="METAL Meta-analysis P value">' >> ${unsorted_Metal_vcf}
echo -e '##INFO=<ID=N_METAL,Number=1,Type=Integer,Description="METAL: Number of samples with marker present">' >> ${unsorted_Metal_vcf}
echo -e '##INFO=<ID=EAF_METAL,Number=A,Type=Float,Description="METAL Average effect allele frequency">' >> ${unsorted_Metal_vcf}
echo -e '##INFO=<ID=BETA_METAL,Number=A,Type=Float,Description="METAL effect estimate">' >> ${unsorted_Metal_vcf}
echo -e '##INFO=<ID=SE_METAL,Number=A,Type=Float,Description="METAL SE of effect estimate">' >> ${unsorted_Metal_vcf}
echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' >> ${unsorted_Metal_vcf}
	# Meta-GWAS summary 
awk 'NR>1 && $1~/^chr[0-9]+/ {print $1"\t"$4"\t"$8"\t"$9"\t"$10"\t"$NF}' ${output_meta}1.tsv | sed 's:\::\t:g' | awk '{printf $1"\t"$2"\t"$1":"$2":"$3":"$4"\t"$3"\t"$4"\t.\t.\tEAF_METAL="$5";BETA_METAL="$6";SE_METAL="$7";P_METAL="$8";N_METAL="$9"\n"}' >> ${unsorted_Metal_vcf} && rm ${output_meta}1.tsv

	# 
bgzip -f ${unsorted_Metal_vcf} 

# no index for unsorted file tabix -f -p vcf ${unsorted_Metal_vcf}.gz

	## sort
zcat ${unsorted_Metal_vcf}.gz | vcf-sort -c > ${sorted_Metal_vcf}

bgzip -f ${sorted_Metal_vcf} && tabix -f -p vcf ${sorted_Metal_vcf}.gz

##########################
######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

