#!/bin/bash

source /broad/software/scripts/useuse

use Tabix

#####
# while read pheno; do echo -e "/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/chr1_22.MultiANC.ukb200k_N193342.has${pheno}.regenie.noNA.tsv\n/broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb/chr1_22.MultiANC.ukb250k_N243350.has${pheno}.regenie.noNA.tsv" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.MultiANC.ukbb450k.ukb200k_N193342.ukb250k_N243350.txt; qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/tmpdir -pe smp 4 -binding linear:4 -l h_vmem=30G -l h_rt=10:00:00 -N MultiANC.ukb200k_250k.${pheno} /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/01.GWAMA.Meta_Analysis.ukb450k.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/has${pheno}.MultiANC.ukbb450k.ukb200k_N193342.ukb250k_N243350.txt /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/chr1_22.has${pheno}.MultiANC.ukbb450k.ukb200k_N193342.ukb250k_N243350; done < <(echo -e "CH\nCHvaf10\nDNMT3A\nTET2\nASXL1\nSF\nDDR")

############### Task: Worked
## 1. Run Meta-Analysis
## 2. sorted summary stats MAF 0.1 and N_studies>1
###########################

######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
##################################
## while read traits
# do
# zcat /medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/chr1_22.has${traits}.21Aug2021_ukbEUR.tsv.gz | awk '(NR==1){print "SNPID\tN\tREF\tALT\tAAF\tBETA\tSE\tP"}(NR>1 && $2~/^[0-9]+/){print $1"\t"($14+$18)"\t"$4"\t"$5"\t"$13"\t"$22"\t"$23"\t"$12}' | awk '!seen[$1]++'| gzip -c > /medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/chr1_22.ukb200k.hg37_eur_${traits}.tsv.gz
	##
# zcat /medpop/esp/mesbah/GWAS_CHIP/inputUKB/ukb450k/250k/hg37/chr1_22.ukb250k_WB.has${traits}.tsv.gz | awk '(NR==1){print "SNPID\tN\tREF\tALT\tAAF\tBETA\tSE\tP"}(NR>1 && $2~/^[0-9]+/){print $1"\t"($14+$18)"\t"$4"\t"$5"\t"$13"\t"$22"\t"$23"\t"$12}' |awk '!seen[$1]++'| gzip -c > /medpop/esp/mesbah/GWAS_CHIP/inputUKB/ukb450k/250k/hg37/chr1_22.ukb250k.hg37_eur_${traits}.regenie.tsv.gz
	## 
# echo -e "/medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/chr1_22.ukb200k.hg37_eur_${traits}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/inputUKB/ukb450k/250k/hg37/chr1_22.ukb250k.hg37_eur_${traits}.regenie.tsv.gz" > /medpop/esp/mesbah/GWAS_CHIP/inputUKB/ukb200k.ukb250k.hg37_eur_${traits}.gwama_input.txt
# done < <(echo -e "CHIP\nDNMT3A\nTET2")
#####################
## EUR UKB 200k, 250k GWAMA for finemapping
## while read trait; do qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/eur_only_summary/ukb450k_eur -pe smp 1 -binding linear:1 -l h_vmem=30G -l h_rt=10:00:00 -N ukb200k_250k.eur_metaGWAS.${trait} /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/01.GWAMA.Meta_Analysis.ukb450k.sh /medpop/esp/mesbah/GWAS_CHIP/inputUKB/ukb200k.ukb250k.hg37_eur_${trait}.gwama_input.txt /broad/hptmp/mesbah/gwas/ukb200k_ukb250k.meta_${trait}_hg37_EUR.gwama /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/eur_only_summary/ukb450k_eur ; done < <(echo -e "CHIP\nDNMT3A\nTET2")
###
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
# awk '(NR==1){print $0}(NR>1 && $15==2){print $0}' ${output_prefix}.out | gzip -c > ${outDir}/$(basename ${output_prefix}.out).gz 

gzip ${output_prefix}.out

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

