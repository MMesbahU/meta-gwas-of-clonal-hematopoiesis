#!/bin/bash


## Note: Run /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/prepGWASsummary.UKB.sh 
# to merge GWAS

### UKB450k WES 
## while read pheno; do qsub -R y -wd /broad/hptmp/mesbah/gwas/ukb450k/450k.step2/tmpdir -N Liftover_allukb450k${pheno} -pe smp 1 -binding linear:1 -l h_rt=20:00:00 -l h_vmem=40G /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/01.PrepSummary/02.prepUKBB_summary.sh /medpop/esp2/mesbah/tools/LiftOVer/picard.jar /broad/hptmp/mesbah/gwas/ukb450k/450k.step2/chr1_22.has${pheno}.regenie.gz /broad/hptmp/mesbah/gwas/ukb450k/450k.step2/tmp.chr1_22.has${pheno}.ukb450k.vcf /broad/hptmp/mesbah/gwas/ukb450k/450k.step2/lifted_hg38.chr1_22.has${pheno}.ukb450k.vcf.gz /medpop/esp2/mesbah/tools/LiftOVer/b37ToHg38.over.chain /broad/hptmp/mesbah/gwas/ukb450k/450k.step2/rejected_lifted_hg38.chr1_22.has${pheno}.ukb450k.vcf.gz /medpop/esp/skoyama/publicdata/reference/Homo_sapiens_assembly38.fasta /broad/hptmp/mesbah/gwas/ukb450k/450k.step2/tmp_info.chr1_22.has${pheno}.ukb450k /broad/hptmp/mesbah/gwas/ukb450k/450k.step2/lifted_hg38.chr1_22.has${pheno}.ukb450k.INFO.tsv; done < <(echo -e "CHIP\nDNMT3A\nTET2\nASXL1")

## UKB 250k GWAS
## while read pheno; do qsub -R y -wd /broad/hptmp/mesbah/gwas/ukb450k/250k.step2/tmpdir -N Liftover_allukb250k${pheno} -pe smp 1 -binding linear:1 -l h_rt=20:00:00 -l h_vmem=40G /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/01.PrepSummary/02.prepUKBB_summary.sh /medpop/esp2/mesbah/tools/LiftOVer/picard.jar /broad/hptmp/mesbah/gwas/ukb450k/250k.step2/chr1_22.ukb250k.has${pheno}.tsv.gz /broad/hptmp/mesbah/gwas/ukb450k/250k.step2/ukb250k_gwama_input/tmp.chr1_22.has${pheno}.ukb250k.vcf /broad/hptmp/mesbah/gwas/ukb450k/250k.step2/ukb250k_gwama_input/lifted_hg38.chr1_22.has${pheno}.ukb250k.vcf.gz /medpop/esp2/mesbah/tools/LiftOVer/b37ToHg38.over.chain /broad/hptmp/mesbah/gwas/ukb450k/250k.step2/ukb250k_gwama_input/rejected_lifted_hg38.chr1_22.has${pheno}.ukb250k.vcf.gz /medpop/esp/skoyama/publicdata/reference/Homo_sapiens_assembly38.fasta /broad/hptmp/mesbah/gwas/ukb450k/250k.step2/ukb250k_gwama_input/tmp_info.chr1_22.has${pheno}.ukb250k /broad/hptmp/mesbah/gwas/ukb450k/250k.step2/ukb250k_gwama_input/lifted_hg38.chr1_22.has${pheno}.ukb250k.INFO.tsv; done < <(echo -e "CHIP\nDNMT3A\nTET2\nASXL1")
## WB UKB 250k GWAS
# while read pheno; do qsub -R y -wd /broad/hptmp/mesbah/gwas/ukb450k/250k.step2/tmpdir -N Liftover_WB_ukb250k.${pheno} -pe smp 1 -binding linear:1 -l h_rt=20:00:00 -l h_vmem=40G /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/01.PrepSummary/02.prepUKBB_summary.sh /medpop/esp2/mesbah/tools/LiftOVer/picard.jar /broad/hptmp/mesbah/gwas/ukb450k/250k.step2/ukb250k_WB/chr1_22.ukb250k_WB.has${pheno}.tsv.gz /broad/hptmp/mesbah/gwas/ukb450k/250k.step2/ukb250k_WB/tmp.chr1_22.has${pheno}.ukb250k_WB.vcf /broad/hptmp/mesbah/gwas/ukb450k/250k.step2/ukb250k_WB/lifted_hg38.chr1_22.has${pheno}.ukb250k_WB.vcf.gz /medpop/esp2/mesbah/tools/LiftOVer/b37ToHg38.over.chain /broad/hptmp/mesbah/gwas/ukb450k/250k.step2/ukb250k_WB/rejected_lifted_hg38.chr1_22.has${pheno}.ukb250k_WB.vcf.gz /medpop/esp/skoyama/publicdata/reference/Homo_sapiens_assembly38.fasta /broad/hptmp/mesbah/gwas/ukb450k/250k.step2/ukb250k_WB/tmp_info.chr1_22.has${pheno}.ukb250k_WB /broad/hptmp/mesbah/gwas/ukb450k/250k.step2/ukb250k_WB/lifted_hg38.chr1_22.has${pheno}.ukb250k_WB.INFO.tsv; done < <(echo -e "CHIP\nDNMT3A\nTET2")

######################################################

############### 200k WES
## All UKB GWAS
# while read pheno; do qsub -R y -wd /broad/hptmp/mesbah/tmpdir -N Liftover_allukb${pheno} -pe smp 1 -binding linear:1 -l h_rt=20:00:00 -l h_vmem=20G /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/01.PrepSummary/02.prepUKBB_summary.sh /medpop/esp2/mesbah/tools/LiftOVer/picard.jar /medpop/esp/mesbah/GWAS_CHIP/inputUKB/allUKB/chr1_22.has${pheno}.11Aug2021_ukb200k.tsv.gz /broad/hptmp/mesbah/ukb_chip/allUKB/tmp.chr1_22.has${pheno}.11Aug2021_ukb200k.vcf /medpop/esp/mesbah/GWAS_CHIP/inputUKB/allUKB/lifted_hg38.chr1_22.has${pheno}.11Aug2021_ukb200k.vcf.gz /medpop/esp2/mesbah/tools/LiftOVer/b37ToHg38.over.chain /medpop/esp/mesbah/GWAS_CHIP/inputUKB/allUKB/rejected_lift_hg38.chr1_22.has${pheno}.11Aug2021_ukb200k.vcf.gz /medpop/esp/skoyama/publicdata/reference/Homo_sapiens_assembly38.fasta /broad/hptmp/mesbah/ukb_chip/allUKB/tmp_info.chr1_22.has${pheno}.11Aug2021_ukb200k /medpop/esp/mesbah/GWAS_CHIP/inputUKB/allUKB/lifted_hg38.chr1_22.has${pheno}.11Aug2021_ukb200k.INFO.tsv; done < <(echo -e "CHIP\nDNMT3A\nTET2")
### EUR only GWAS
## while read pheno; do qsub -R y -wd /broad/hptmp/mesbah/tmpdir -N Liftover_ukb_eur${pheno} -pe smp 1 -binding linear:1 -l h_rt=20:00:00 -l h_vmem=20G /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/01.PrepSummary/02.prepUKBB_summary.sh /medpop/esp2/mesbah/tools/LiftOVer/picard.jar /medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/chr1_22.has${pheno}.21Aug2021_ukbEUR.tsv.gz /broad/hptmp/mesbah/ukb_chip/eurUKB/tmp.chr1_22.has${pheno}.21Aug2021_ukbEUR.vcf /medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/lifted_hg38.chr1_22.has${pheno}.21Aug2021_ukbEUR.vcf.gz /medpop/esp2/mesbah/tools/LiftOVer/b37ToHg38.over.chain /medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/rejected_lift_hg38.chr1_22.has${pheno}.21Aug2021_ukbEUR.vcf.gz /medpop/esp/skoyama/publicdata/reference/Homo_sapiens_assembly38.fasta /broad/hptmp/mesbah/ukb_chip/eurUKB/tmp_info.chr1_22.has${pheno}.21Aug2021_ukbEUR /medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/lifted_hg38.chr1_22.has${pheno}.21Aug2021_ukbEUR.INFO.tsv; done < <(echo -e "CHIP\nDNMT3A\nTET2")
################

source /broad/software/scripts/useuse

use Picard-Tools
use Tabix
use VCFtools

######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
#########

### input files
picard=${1}
inGWAS=${2}
tmpVCF=${3}
outVCF=${4}
chain=${5}
rejectVCF=${6}
refTarget=${7}
out_info_file=${8}
final_output=${9}


## Prepare input VCF

# zcat /medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/chr1_22.hasCHIP.21Aug2021_ukbEUR.tsv.gz | head -2
# Name	Chr	Pos_hg37	Ref	Alt	Trait	Cohort	Model	Effect	LCI_Effect UCI_Effect	Pval	AAF	Num_Cases	Cases_Ref	Cases_Het	Cases_Alt	Num_Controls	Controls_Ref	Controls_Het	Controls_Alt	REGENIE_BETA	REGENIE_SE INFO	MAC
# rs367896724	1	10177	A	AC	hasCHIP	eur_ukb	ADD-WGR-FIRTH	0.990115	0.947629	1.03451	0.657099	0.398959	9789	3054	5965	770	157924	48698	96715	12511	-0.00993377	0.0223773	0.467961	133822

echo -e '##fileformat=VCFv4.2\n##INFO=<ID=PVAL,Number=A,Type=Float,Description="GWAS Pvalue">\n##INFO=<ID=AAF,Number=A,Type=Float,Description="Alternate Allele Frequency">\n##INFO=<ID=NCASES,Number=1,Type=Integer,Description="Number of CHIP Cases">\n##INFO=<ID=NCONTROL,Number=1,Type=Integer,Description="Number of CHIP Controls">\n##INFO=<ID=BETA,Number=A,Type=Float,Description="SNP effect estimate beta">\n##INFO=<ID=SE,Number=A,Type=Float,Description="SE estimate of SNP effect">\n##INFO=<ID=INFO,Number=A,Type=Float,Description="Imputation r square">\n##INFO=<ID=MAC,Number=1,Type=Integer,Description="Minor allele count">\n##INFO=<ID=N,Number=A,Type=Float,Description="Total Samples with genotype">\n##INFO=<ID=old_hg19_id,Number=1,Type=String,Description="UKB HG19/hg37 coordinates">\n##INFO=<ID=old_ukbrsid,Number=1,Type=String,Description="UKB rsids">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' > ${tmpVCF}

## Keep Chr1-22 and numeric positions
zcat ${inGWAS} | awk 'NR>1 && $2~/[0-9]+/ && $3~/[0-9]+/ {print $2"\t"$3"\t"$1"\t"$4"\t"$5"\t.\t.\tPVAL="$12";AAF="$13";BETA="$22";SE="$23";INFO="$24";MAC="$25";NCASES="$14";NCONTROL="$18";old_hg19_id="$2":"$3":"$4":"$5";old_ukbrsid="$1";N="($14+$18)}' >> ${tmpVCF}

## Run Picard to liftover
java -Xmx40G -jar ${picard} LiftoverVcf \
	-I ${tmpVCF} \
	-O ${outVCF} \
	-C ${chain} \
	--MAX_RECORDS_IN_RAM 100000 \
	--REJECT ${rejectVCF} \
	-R ${refTarget}

## Prepare Output summary File for meta-analysis
# rm ${tmpVCF}

## Keep chr1 to chr 22
# $(echo chr{1..22} chrX)
## --gzvcf ${outVCF}
tabix -h ${outVCF} $(echo chr{1..22}) | vcftools \
	--vcf - \
	--get-INFO AAF \
	--get-INFO BETA \
	--get-INFO SE \
	--get-INFO PVAL \
	--get-INFO N \
	--get-INFO old_hg19_id \
	--out ${out_info_file}

## Add CHROM:POS:REF:ALT 
head -1 ${out_info_file}.INFO | awk '{print "SNPID\t"$0}' > ${final_output}

# Exclude Duplicate SNPIDs
awk 'NR>1{print $1":"$2":"$3":"$4"\t"$0}' ${out_info_file}.INFO |  awk '!seen[$1]++' | sort -k1 -V >> ${final_output}

rm ${out_info_file}.INFO

# compress 
bgzip -f ${final_output}

tabix -f -s 2 -b 3 -e 3 ${final_output}.gz

#####

######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########


