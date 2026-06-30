#!/bin/bash

## Sept 06
# while read lines; do qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/SumHer -N SumHer.Sept6_pop$(echo ${lines} | awk '{print $2}' )_sample$(echo ${lines} | awk '{print $3}' ).eur_metaGWAS.$(echo ${lines} | awk '{print $1}' ) -l h_rt=20:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.SNP_heritability/01_1.run_ldak.sh /broad/hptmp/mesbah/gwas/sumher/ldak.lifted_hg37.eur_metaGWAS.$(echo ${lines} | awk '{print $1}' ).GWAMA.hg37_dbSNP.txt /broad/hptmp/mesbah/gwas/sumher/ref/bld_ldak/bld.ldak.genotyped.gbr.tagging /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/SumHer/sumher.pop$(echo ${lines} | awk '{print $2}' )_sample$(echo ${lines} | awk '{print $3}' ).eur_metaGWAS.$(echo ${lines} | awk '{print $1}' ) $(echo ${lines} | awk '{print $2}' ) $(echo ${lines} | awk '{print $3}' ); done < <(echo -e "CHIP\t0.0631\t0.0639\nDNMT3A\t0.0359\t0.0349\nTET2\t0.0131\t0.0137" )
# 
### multi-ancestry:
# while read lines; do qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/SumHer -N SumHer.Sept6_pop$(echo ${lines} | awk '{print $2}' )_sample$(echo ${lines} | awk '{print $3}' ).multiAncestry_metaGWAS.$(echo ${lines} | awk '{print $1}' ) -l h_rt=20:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.SNP_heritability/01_1.run_ldak.sh /broad/hptmp/mesbah/gwas/sumher/ldak.lifted_hg37.metaGWAS.$(echo ${lines} | awk '{print $1}' ).TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.hg37_dbSNP.txt /broad/hptmp/mesbah/gwas/sumher/ref/bld_ldak/bld.ldak.genotyped.gbr.tagging /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/SumHer/sumher.pop$(echo ${lines} | awk '{print $2}' )_sample$(echo ${lines} | awk '{print $3}' ).multiAncestry_metaGWAS.$(echo ${lines} | awk '{print $1}' ) $(echo ${lines} | awk '{print $2}' ) $(echo ${lines} | awk '{print $3}' ); done < <(echo -e "CHIP\t0.0625\t0.0631\nDNMT3A\t0.0355\t0.0331\nTET2\t0.013\t0.0129" )

######

### Aug 25, 2022
# EUR GWAS
# for files in $(ls -l /broad/hptmp/mesbah/gwas/sumher/ldak.lifted_hg37.eur_metaGWAS.*.GWAMA.hg37_dbSNP.txt | awk '{print $NF}'); do qsub -R y -wd /broad/hptmp/mesbah/gwas/sumher -N eur_SumHer.aug25 -l h_rt=20:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.SNP_heritability/01_1.run_ldak.sh ${files} /broad/hptmp/mesbah/gwas/sumher/ref/bld_ldak/bld.ldak.genotyped.gbr.tagging /broad/hptmp/mesbah/gwas/sumher/sumher.pop05_eur05.$(basename ${files} ".txt") 0.05 0.05; done

# Multi-ancestry GWAS
# for files in $(ls -l ldak.lifted_hg37.metaGWAS.*.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.hg37_dbSNP.txt | awk '{print $NF}'); do qsub -R y -wd /broad/hptmp/mesbah/gwas/sumher -N multiAncestry_SumHer.aug25 -l h_rt=20:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.SNP_heritability/01_1.run_ldak.sh ${files} /broad/hptmp/mesbah/gwas/sumher/ref/bld_ldak/bld.ldak.genotyped.gbr.tagging /broad/hptmp/mesbah/gwas/sumher/sumher.pop05_multiAncestry05.$(basename ${files} ".txt") 0.05 0.05; done
##########

### March 8, 2022
## mixed GWAS: for files in $(ls -l /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/in_ldak_eaf001_min2Studies.8March2022.lifted_hg37.GWAMA.mixedSamples.meta_*.txt | awk '{print $NF}'); do qsub -R y -wd /broad/hptmp/mesbah/ukb_chip/tmpdir -N meta_SumHer.8March2022 -l h_rt=20:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/09.Heritability/1.run_ldak.sh ${files} /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/bld.ldak.genotyped.gbr.tagging /broad/hptmp/mesbah/ukb_chip/h2/output_ldak/sumher.pop05_meta05.$(basename ${files} ".txt") 0.05 0.05; done
## Meta4: files=/broad/hptmp/mesbah/ukb_chip/h2/input_ldak/in_ldak_eaf001_min2Studies.8March2022.lifted_hg37.GWAMA.meta_CHIP.allsamples_topmed_ukbb_mgbb_BioVU.txt; qsub -R y -wd /broad/hptmp/mesbah/ukb_chip/tmpdir -N chip_meta4.SumHer.8March2022 -l h_rt=20:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/09.Heritability/1.run_ldak.sh ${files} /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/bld.ldak.genotyped.gbr.tagging /broad/hptmp/mesbah/ukb_chip/h2/output_ldak/sumher.pop05_meta05.$(basename ${files} ".txt") 0.05 0.05 

## for files in $(ls -lhrt /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/in_ldak_eaf001_min2Studies.* | awk '{print $NF}' | tail -3); do qsub -R y -wd /broad/hptmp/mesbah/ukb_chip/tmpdir -N meta_SumHer.8March2022 -l h_rt=20:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/09.Heritability/1.run_ldak.sh ${files} /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/bld.ldak.genotyped.gbr.tagging /broad/hptmp/mesbah/ukb_chip/h2/output_ldak/sumher.pop05_meta05.$(basename ${files} ".txt") 0.05 0.05; done

###

## for files in $(ls -l /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/summary_for_ldak.chr1_22.has*.21Aug2021_ukbEUR.tsv | awk '{print $NF}'); do qsub -R y -wd /broad/hptmp/mesbah/ukb_chip/tmpdir -N SumHer -l h_rt=20:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/09.Heritability/1.run_ldak.sh ${files} /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/bld.ldak.genotyped.gbr.tagging /broad/hptmp/mesbah/ukb_chip/h2/output_ldak/sumher.allukb059_eur052.$(basename ${files} ".tsv") 0.059 0.052; done

## EUR TOPMED
## qsub -R y -wd /broad/hptmp/mesbah/ukb_chip/tmpdir -N SumHer.topmedEUR -l h_rt=20:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/09.Heritability/1.run_ldak.sh /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/summary_for_ldak.topmed.hg37.CHIP.EUR.results.txt /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/bld.ldak.genotyped.gbr.tagging /broad/hptmp/mesbah/ukb_chip/h2/output_ldak/sumher.alltopmed10_eur10.summary_for_ldak.topmed.hg37.CHIP.EUR.results 0.1 0.1

## EUR 4 studies
# qsub -R y -wd /broad/hptmp/mesbah/ukb_chip/tmpdir -N SumHer.eurmeta4 -l h_rt=20:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/09.Heritability/1.run_ldak.sh /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/summary_for_ldak.lifted_hg37.GWAMA.EUR_only_CHIP_GWAS.TopMedUKbbMGbbBioVU.eaf001_min2Studies.txt /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/bld.ldak.genotyped.gbr.tagging /broad/hptmp/mesbah/ukb_chip/h2/output_ldak/sumher.allmeta410_eur10.lifted_hg37.GWAMA.EUR_only_CHIP_GWAS.TopMedUKbbMGbbBioVU.eaf001_min2Studies 0.1 0.1
# qsub -R y -wd /broad/hptmp/mesbah/ukb_chip/tmpdir -N SumHer.eurmeta4 -l h_rt=20:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/09.Heritability/1.run_ldak.sh /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/summary_for_ldak.lifted_hg37.GWAMA.EUR_only_CHIP_GWAS.TopMedUKbbMGbbBioVU.eaf001_min2Studies.txt /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/bld.ldak.genotyped.gbr.tagging /broad/hptmp/mesbah/ukb_chip/h2/output_ldak/sumher.allmeta4pop05_eur05.lifted_hg37.GWAMA.EUR_only_CHIP_GWAS.TopMedUKbbMGbbBioVU.eaf001_min2Studies 0.05 0.05


####### HapMap tagging
# qsub -R y -wd /broad/hptmp/mesbah/ukb_chip/tmpdir -N SumHer_HapMapTag.eurmeta4 -l h_rt=20:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/09.Heritability/1.run_ldak.sh /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/summary_for_ldak.lifted_hg37.GWAMA.EUR_only_CHIP_GWAS.TopMedUKbbMGbbBioVU.eaf001_min2Studies.txt /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/bld.ldak.hapmap.gbr.tagging /broad/hptmp/mesbah/ukb_chip/h2/output_ldak/sumher_HapMapTag.allmeta4pop05_eur05.lifted_hg37.GWAMA.EUR_only_CHIP_GWAS.TopMedUKbbMGbbBioVU.eaf001_min2Studies 0.05 0.05

# qsub -R y -wd /broad/hptmp/mesbah/ukb_chip/tmpdir -N SumHer_HapMapTag.topmed -l h_rt=20:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/09.Heritability/1.run_ldak.sh /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/summary_for_ldak.topmed.hg37.CHIP.EUR.results.txt /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/bld.ldak.hapmap.gbr.tagging /broad/hptmp/mesbah/ukb_chip/h2/output_ldak/sumher_HapMapTag.allTopmedpop05_eur05.summary_for_ldak.topmed.hg37.CHIP.EUR.results 0.05 0.05

# qsub -R y -wd /broad/hptmp/mesbah/ukb_chip/tmpdir -N SumHer_HapMapTag.ukb -l h_rt=20:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/09.Heritability/1.run_ldak.sh /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/summary_for_ldak.chr1_22.hasCHIP.21Aug2021_ukbEUR.tsv /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/bld.ldak.hapmap.gbr.tagging /broad/hptmp/mesbah/ukb_chip/h2/output_ldak/sumher_HapMapTag.allUKBpop05_eur05.chr1_22.hasCHIP.21Aug2021_ukbEUR 0.05 0.05
##############

## Get the software
# wget https://dougspeed.com/wp-content/uploads/ldak5.2.linux_.zip
# unzip ldak5.2.linux_.zip
## HapMap tagging: bld.ldak.hapmap.gbr.tagging

#### http://dougspeed.com/pre-computed-tagging-files/
# Get precomputed Tagging SNPs
# LDAK-Thin Model: wget https://genetics.ghpc.au.dk/doug/ldak.thin.genotyped.gbr.tagging.gz
# BLD-LDAK-Lite-Alpha Model: wget https://genetics.ghpc.au.dk/doug/bld.ldak.lite.alpha.hapmap.gbr.tagging.gz
# BLD-LDAK Model: wget -P /broad/hptmp/mesbah/ukb_chip/ https://genetics.ghpc.au.dk/doug/bld.ldak.genotyped.gbr.tagging.gz


## Sumher: http://dougspeed.com/snp-heritability/
summary_LDAK_format=${1} # to specify the file containing the summary statistics.
taggFile=${2}  # to specify the tagging file
out_file_prefix=${3}
pop_prev=${4} # prevalence in general population
gwas_prev=${5} # prevalence in GWAS population

### If the summary statistics come from analysing a binary phenotype, then you can use --prevalence <float> and --ascertainment <float> to specify the proportion of cases in the population and in the GWAS; LDAK will then also report estimates of variance explained on the liability scale.

######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
#########
/medpop/esp2/mesbah/tools/LDAK5.2/ldak5.2.linux \
	--summary ${summary_LDAK_format} \
	--tagfile ${taggFile} \
	--sum-hers ${out_file_prefix} \
	--prevalence ${pop_prev} \
	--ascertainment ${gwas_prev} \
	--check-sums NO

######################

######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########


