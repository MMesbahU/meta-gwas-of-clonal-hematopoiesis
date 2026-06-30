#!/bin/bash

## EUR: while read trait; do qsub -R y -wd /broad/hptmp/mesbah/gwas/ukb450k/metagwas/tmpdir -N annot_gwama_EUR.${trait} -pe smp 1 -binding linear:1 -l h_rt=20:00:00 -l h_vmem=30G /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/03.GWAMA_VCF.sh /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/eur_metaGWAS.${trait}.*_mgb_biovu.out.gz /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/tmp.eur_metaGWAS.${trait}.GWAMA.vcf /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/sorted.eur_metaGWAS.${trait}.GWAMA.vcf /broad/hptmp/mesbah/gwas/ukb450k/metagwas/sorted_metal_eur.topmed_ukb450k_mgbb53k_biovu54k.${trait}.vcf.gz /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/eur_metaGWAS.${trait}.GWAMA.hg38.vcf.gz 'P_METAL,N_METAL,EAF_METAL,BETA_METAL,SE_METAL' /broad/hptmp/mesbah/RefSeq/dbSNP/hg38/chr1_X.hg38.GCF_000001405.39.vcf.gz /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/eur_metaGWAS.${trait}.GWAMA.hg38_dbSNP.vcf.gz; done < <(echo -e "CHIP\nDNMT3A\nTET2")

# Run: while read trait; do qsub -R y -wd /broad/hptmp/mesbah/gwas/ukb450k/metagwas/tmpdir -N annot_gwama_${trait} -pe smp 1 -binding linear:1 -l h_rt=20:00:00 -l h_vmem=30G /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/03.GWAMA_VCF.sh /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/metaGWAS.${trait}.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.out.gz /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/tmp.metaGWAS.${trait}.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.vcf /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/sorted.metaGWAS.${trait}.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.vcf /broad/hptmp/mesbah/gwas/ukb450k/metagwas/sorted_metal.topmed64k_ukb450k_mgbb53k_biovu54k.${trait}.vcf.gz /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/metaGWAS.${trait}.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.hg38.vcf.gz 'P_METAL,N_METAL,EAF_METAL,BETA_METAL,SE_METAL' /broad/hptmp/mesbah/RefSeq/dbSNP/hg38/chr1_X.hg38.GCF_000001405.39.vcf.gz /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/metaGWAS.${trait}.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.hg38_dbSNP.vcf.gz; done < <(echo -e "CHIP\nDNMT3A\nTET2")
###
source /broad/software/scripts/useuse

use Picard-Tools

use Tabix

use Bcftools

use VCFtools
 
######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
#########
gwama_file=${1} # gziped GWAMA file

tmpVCF=${2} # tmp VCF 

sorted_tmpVCF=${3} # sorted tmp VCF

METAL_VCF=${4} # METAL annotation VCF

annot_VCF=${5} # METAL "N" annotated GWAMA VCF

annotations=${6} # 'P_METAL,N_METAL,EAF_METAL,BETA_METAL,SE_METAL' 

dbSNP_hg38=${7}

outVCFhg38=${8}

# outVCFhg37=${5}
# tmp_outVCFhg37=${6}

# dbSNP_hg38=${7}
# dbSNP_hg37=${8}

# hg38tohg19=${9}
# Picard=${10}
# rejectVCF=${11}
# refTarget=${12}

## convert to vcf
## https://genomics.ut.ee/en/tools/gwama-tutorial
echo -e '##fileformat=VCFv4.2' > ${tmpVCF}
echo -e '##INFO=<ID=varID_hg38,Number=1,Type=String,Description="hg38 chr:pos:ref:alt ids">' >> ${tmpVCF}
echo -e '##INFO=<ID=EffectAllele,Number=1,Type=String,Description="Effect allele (GWAMA reference_allele)">' >> ${tmpVCF}
echo -e '##INFO=<ID=OtherAllele,Number=1,Type=String,Description="GWAMA other allele (non reference allele)">' >> ${tmpVCF}
echo -e '##INFO=<ID=EAF,Number=A,Type=Float,Description="Average effect allele frequency (i.e. freq for reference_allele)">' >> ${tmpVCF}
echo -e '##INFO=<ID=BETA,Number=A,Type=Float,Description="GWAMA effect estimate">' >> ${tmpVCF}
echo -e '##INFO=<ID=SE,Number=A,Type=Float,Description="GWAMA se of effect estimate">' >> ${tmpVCF}
echo -e '##INFO=<ID=BETA_95L,Number=A,Type=Float,Description="Lower 95% CI for BETA">' >> ${tmpVCF}
echo -e '##INFO=<ID=BETA_95U,Number=A,Type=Float,Description="Upper 95% CI for BETA">' >> ${tmpVCF}
echo -e '##INFO=<ID=Z,Number=A,Type=Float,Description="Z-score">' >> ${tmpVCF}
echo -e '##INFO=<ID=Pvalue,Number=A,Type=Float,Description="GWAMA Meta-analysis P value">' >> ${tmpVCF}
echo -e '##INFO=<ID=log10_Pvalue,Number=A,Type=Float,Description="Absolute value of logarithm of meta-analysis p value to the base of 10">' >> ${tmpVCF}
echo -e '##INFO=<ID=Q_stat,Number=A,Type=Float,Description="Cochran heterogeneity statistic">' >> ${tmpVCF} 
echo -e '##INFO=<ID=Q_Pvalue,Number=A,Type=Float,Description="Cochran heterogeneity statistic P value">' >> ${tmpVCF}
echo -e '##INFO=<ID=I_sqr,Number=A,Type=Float,Description="Heterogeneity index I2 by Higgins et al 2003">' >> ${tmpVCF}
echo -e '##INFO=<ID=N_Studies,Number=1,Type=Integer,Description="Number of studies with marker present">' >> ${tmpVCF}
echo -e '##INFO=<ID=N_Samples,Number=1,Type=Integer,Description="Number of samples with marker present (will be NA if marker is present in any input file where N column is not present)">' >> ${tmpVCF}
echo -e '##INFO=<ID=Effect_Direction,Number=1,Type=String,Description="Summary of effect directions in all cohort">' >> ${tmpVCF}
echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' >> ${tmpVCF}

## Fill the vcf 
# zcat  /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/metaGWAS.TET2.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.out.gz | head -2
	# GWAMA Summary
zcat ${gwama_file} | awk 'NR>1 && $1~/^chr[0-9]+/ {print $0}' | sed 's:\::\t:g' | awk '{print $1"\t"$2"\t"$1":"$2":"$3":"$4"\t"$3"\t"$4"\t.\t.\tvarID_hg38="$1":"$2":"$3":"$4";EffectAllele="$5";OtherAllele="$6";EAF="$7";BETA="$8";SE="$9";BETA_95L="$10";BETA_95U="$11";Z="$12";Pvalue="$13";log10_Pvalue="$14";Q_stat="$15";Q_Pvalue="$16";I_sqr="$17";N_Studies="$18";N_Samples="$19";Effect_Direction="$20}' >> ${tmpVCF}

# Sort by chromosome
cat ${tmpVCF} | vcf-sort -c > ${sorted_tmpVCF}

rm ${tmpVCF}

# compress  
bgzip -f ${sorted_tmpVCF}

tabix -f -p vcf ${sorted_tmpVCF}.gz

### Add Sample Size N and P value from METAL 
bcftools annotate \
	-a ${METAL_VCF} \
	-c ${annotations} ${sorted_tmpVCF}.gz \
	-o ${annot_VCF} \
	-Oz
# 
tabix -f -p vcf ${annot_VCF}

#### Add rsids from dbSNP 155
# bcftools annotate -a /broad/hptmp/mesbah/RefSeq/dbSNP/hg38/chr1_X.hg38.GCF_000001405.39.vcf.gz -c ID /medpop/esp/mesbah/GWAMA_VCF/hg38/metaGWAS.CHIP.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.hg38.vcf.gz -Oz -o /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/metaGWAS.CHIP.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.hg38_dbSNP.vcf.gz &
# bcftools annotate \
#	-a ${dbSNP_hg38} \
#	-c ID,GENEINFO,NSF,NSM,NSN,SYN,ASS,DSS,dbSNPBuildID,SSR ${annot_VCF} \
#	-Oz \
#	-o ${outVCFhg38}

# tabix -p vcf ${outVCFhg38}

## Lift to hg37 and add rdids
## Run Picard to liftover
# java -jar ${Picard} LiftoverVcf \
#	-I ${sorted_tmpVCF}.gz \
#	-O ${tmp_outVCFhg37} \
#	-C ${hg38tohg19} \
#	--REJECT ${rejectVCF} \
#	-R ${refTarget}

# tabix -p vcf ${tmp_outVCFhg37}

	##
# bcftools annotate \
#	-a ${dbSNP_hg37} \
#	-c ID ${tmp_outVCFhg37} \
#	-Oz \
#	-o ${outVCFhg37}

# tabix -p vcf ${outVCFhg37}


######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

