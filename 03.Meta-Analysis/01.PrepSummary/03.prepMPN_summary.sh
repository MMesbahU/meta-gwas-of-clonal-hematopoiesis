#!/bin/bash

## qsub -R y -wd /broad/hptmp/mesbah/tmpdir -N Liftover_mpn -pe smp 1 -binding linear:1 -l h_rt=20:00:00 -l h_vmem=20G /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/01.PrepSummary/03.prepMPN_summary.sh /medpop/esp2/mesbah/tools/LiftOVer/picard.jar /medpop/esp/mesbah/GWAS_CHIP/MPN_summary/MPN_metaGWAS_sumstats.tsv.gz /broad/hptmp/mesbah/ukb_chip/MPN/tmp.MPN_metaGWAS_sumstats.vcf /medpop/esp/mesbah/GWAS_CHIP/MPN_summary/lifted_hg38.MPN_metaGWAS_sumstats.vcf.gz /medpop/esp2/mesbah/tools/LiftOVer/b37ToHg38.over.chain /medpop/esp/mesbah/GWAS_CHIP/MPN_summary/rejected_lift_hg38.MPN_metaGWAS_sumstats.vcf.gz /medpop/esp/skoyama/publicdata/reference/Homo_sapiens_assembly38.fasta /broad/hptmp/mesbah/ukb_chip/MPN/tmp_info.MPN_metaGWAS_sumstats /medpop/esp/mesbah/GWAS_CHIP/MPN_summary/lifted_hg38.MPN_metaGWAS_sumstats.INFO.tsv

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
# zcat /medpop/esp/mesbah/GWAS_CHIP/MPN_summary/MPN_metaGWAS_sumstats.tsv.gz | head -2
# MarkerName	RSID	CHR	POS	REF	ALT	Effect	StdErr	pvalue	MAF
# 1:768448_G_A	rs12562034	1	768448	G	A	-0.0253	0.0615	0.681	0.1076
echo -e '##fileformat=VCFv4.2\n##INFO=<ID=pvalue,Number=A,Type=Float,Description="GWAS Pvalue">\n##INFO=<ID=MAF,Number=A,Type=Float,Description="Minor Allele Frequency">\n##INFO=<ID=Effect,Number=A,Type=Float,Description="SNP effect estimate beta">\n##INFO=<ID=StdErr,Number=A,Type=Float,Description="SE estimate of SNP effect">\n##INFO=<ID=MarkerName,Number=1,Type=String,Description="Original MPN hg19/hg37 coordinates">\n##INFO=<ID=RSID,Number=1,Type=String,Description="Old MPN RSID">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' > ${tmpVCF}

# zcat /medpop/esp/mesbah/GWAS_CHIP/MPN_summary/MPN_metaGWAS_sumstats.tsv.gz | head -2
# MarkerName    RSID    CHR     POS     REF     ALT     Effect  StdErr  pvalue  MAF
# 1:768448_G_A  rs12562034      1       768448  G       A       -0.0253 0.0615  0.681   0.1076

zcat ${inGWAS} | awk 'NR>1{print $3"\t"$4"\t"$2"\t"$5"\t"$6"\t.\t.\tpvalue="$9";MAF="$10";Effect="$7";StdErr="$8";MarkerName="$1";RSID="$2}' >> ${tmpVCF}

## Run Picard to liftover
java -jar ${picard} LiftoverVcf \
	-I ${tmpVCF} \
	-O ${outVCF} \
	-C ${chain} \
	--REJECT ${rejectVCF} \
	-R ${refTarget}

## Prepare Output summary File for meta-analysis
rm ${tmpVCF}
##
vcftools \
	--gzvcf ${outVCF} \
	--get-INFO MAF \
	--get-INFO Effect \
	--get-INFO StdErr \
	--get-INFO pvalue \
	--get-INFO MarkerName \
	--get-INFO RSID \
	--out ${out_info_file}

## Add CHROM:POS:REF:ALT 
head -1 ${out_info_file}.INFO | awk '{print "SNPID\t"$0}' > ${final_output}

awk 'NR>1{print $1":"$2":"$3":"$4"\t"$0}' ${out_info_file}.INFO >> ${final_output}

rm ${out_info_file}.INFO

# compress 
bgzip ${final_output}

tabix -s 2 -b 3 -e 3 ${final_output}.gz

#####

######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

