#!/bin/bash


### Tasks:
# 1. Prepare VCF
# 2. liftover hg38 to hg37
##

## qsub -R y -wd /broad/hptmp/mesbah/dataset/ch_gwas/tmpdir -t 1-$(wc -l /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/list_of_113_gwama_output_files.list |awk '{print $1}') -pe smp 1 -binding linear:1 -l h_vmem=30G -l h_rt=10:00:00 -N hg38_to_hg37 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/11.etcs/GWAMA_hg38_to_hg37.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/list_of_113_gwama_output_files.list /medpop/esp2/mesbah/tools/LiftOVer/grch38_to_grch37.over.chain /medpop/esp2/mesbah/tools/LiftOVer/picard.jar /broad/hptmp/mesbah/RefSeq/human_g1k_v37.fasta /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/out/hg38 /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/out/hg37

#### Note: ## SORT step in LiftoverVcf requires lots of Memory!!!! 
# solution: DISABLE_SORT=true
#####

### Prepare Reference File
###### wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.{fasta.gz,fasta.fai}
## Get reference from MOCHA
# wget https://software.broadinstitute.org/software/mocha/mocha.GRCh37.zip
# unzip mocha.GRCh37.zip human_g1k_v37.{fasta,fasta.fai}
# prepare Dict file
# use Picard-Tools; java -jar /medpop/esp2/mesbah/tools/LiftOVer/picard.jar CreateSequenceDictionary R=/broad/hptmp/mesbah/RefSeq/human_g1k_v37.fasta O=/broad/hptmp/mesbah/RefSeq/human_g1k_v37.dict
###################################
### get list of meta summary files:
## ls  /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/*.out.gz | awk '{print $NF}' | sort -V  > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/list_of_113_gwama_output_files.list
## 

###############################
source /broad/software/scripts/useuse

use Picard-Tools

use Tabix

use Bcftools

use VCFtools
######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
#########
GWAMA_file_list=${1}

hg38tohg37=${2} # /medpop/esp2/mesbah/tools/LiftOVer/grch38_to_grch37.over.chain

Picard=${3} # /medpop/esp2/mesbah/tools/LiftOVer/picard.jar

refTarget=${4} #/broad/hptmp/mesbah/RefSeq/human_g1k_v37.fasta

hg38VCFDir=${5}

liftedOutDir=${6}
###
## input files
gwama_file=$(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${GWAMA_file_list} )

inVCF=${hg38VCFDir}/hg38.$(basename $(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${GWAMA_file_list} ) ".out.gz").vcf

sorted_VCF=${hg38VCFDir}/hg38.$(basename $(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${GWAMA_file_list} ) ".out.gz").sorted.vcf

##
liftedVCF=${liftedOutDir}/hg37.$(basename $(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${GWAMA_file_list} ) ".out.gz").vcf.gz

rejectVCF=${liftedOutDir}/rejected_hg37.$(basename $(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${GWAMA_file_list} ) ".out.gz").vcf.gz

#############################
### Prepare VCF file
echo -e "input file: ${gwama_file} \n"

## convert to vcf
## https://genomics.ut.ee/en/tools/gwama-tutorial
echo -e '##fileformat=VCFv4.2' > ${inVCF}
echo -e '##INFO=<ID=varID_hg38,Number=1,Type=String,Description="hg38 chr:pos:ref:alt ids">' >> ${inVCF}
echo -e '##INFO=<ID=EffectAllele,Number=1,Type=String,Description="Effect allele (GWAMA reference_allele)">' >> ${inVCF}
echo -e '##INFO=<ID=OtherAllele,Number=1,Type=String,Description="GWAMA other allele (non reference allele)">' >> ${inVCF}
echo -e '##INFO=<ID=EAF,Number=A,Type=Float,Description="Average effect allele frequency (i.e. freq for reference_allele)">' >> ${inVCF}
echo -e '##INFO=<ID=BETA,Number=A,Type=Float,Description="GWAMA effect estimate">' >> ${inVCF}
echo -e '##INFO=<ID=SE,Number=A,Type=Float,Description="GWAMA se of effect estimate">' >> ${inVCF}
echo -e '##INFO=<ID=BETA_95L,Number=A,Type=Float,Description="Lower 95% CI for BETA">' >> ${inVCF}
echo -e '##INFO=<ID=BETA_95U,Number=A,Type=Float,Description="Upper 95% CI for BETA">' >> ${inVCF}
echo -e '##INFO=<ID=Z,Number=A,Type=Float,Description="Z-score">' >> ${inVCF}
echo -e '##INFO=<ID=Pvalue,Number=A,Type=String,Description="GWAMA Meta-analysis P value">' >> ${inVCF}
echo -e '##INFO=<ID=log10_Pvalue,Number=A,Type=Float,Description="Absolute value of logarithm of meta-analysis p value to the base of 10">' >> ${inVCF}
echo -e '##INFO=<ID=Q_stat,Number=A,Type=Float,Description="Cochran heterogeneity statistic">' >> ${inVCF} 
echo -e '##INFO=<ID=Q_Pvalue,Number=A,Type=Float,Description="Cochran heterogeneity statistic P value">' >> ${inVCF}
echo -e '##INFO=<ID=I_sqr,Number=A,Type=Float,Description="Heterogeneity index I2 by Higgins et al 2003">' >> ${inVCF}
echo -e '##INFO=<ID=N_Studies,Number=1,Type=Integer,Description="Number of studies with marker present">' >> ${inVCF}
echo -e '##INFO=<ID=N_Samples,Number=1,Type=Integer,Description="Number of samples with marker present (will be NA if marker is present in any input file where N column is not present)">' >> ${inVCF}
echo -e '##INFO=<ID=Effect_Direction,Number=1,Type=String,Description="Summary of effect directions in all cohort">' >> ${inVCF}
echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' >> ${inVCF}

# rs_number	reference_allele	other_allele	eaf	beta	se	beta_95L	beta_95U	p-value	_-log10_p-value	q_statistic	q_p-value	i2	n_studies	n_samples	effects
# chr1:10612:A:C	C	A	0.000051	1.574385	0.756204	0.092225	3.056545	2.081958	0.037366	1.427524	0.000000	1.000000	-nan	1	435506	+????

zcat ${gwama_file} | awk 'BEGIN {OFS = "\t"} NR > 1 && $1 ~ /^chr[0-9]+/ {split($1, fields, ":"); chr = fields[1]; pos = fields[2]; ref = fields[3]; alt = fields[4]; varID = $1; print chr, pos, varID, ref, alt, ".", ".", "varID_hg38=" varID ";EffectAllele=" $2 ";OtherAllele=" $3 ";EAF=" $4 ";BETA=" $5 ";SE=" $6 ";BETA_95L=" $7 ";BETA_95U=" $8 ";Z=" $9 ";Pvalue=" $10 ";log10_Pvalue=" $11 ";Q_stat=" $12 ";Q_Pvalue=" $13 ";I_sqr=" $14 ";N_Studies=" $15 ";N_Samples=" $16 ";Effect_Direction=" $17 }' >> ${inVCF}

## 
# Sort by chromosome
cat ${inVCF} | vcf-sort -c > ${sorted_VCF}

rm ${inVCF}

# compress  
bgzip -f ${sorted_VCF}


############################

#############################
# Run picard LiftoverVcf
# java -Xmx20G -jar /medpop/esp2/mesbah/tools/LiftOVer/picard.jar LiftoverVcf I=/broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/metaGWAS.CHIP.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.hg38.vcf.gz O=/broad/hptmp/mesbah/gwas/lifted_hg37/v2.lifted_hg37.GWAMA_CHIP.vcf.gz REJECT=/broad/hptmp/mesbah/gwas/lifted_hg37/v2.rejected_hg37.GWAMA_CHIP.vcf.gz CHAIN=/medpop/esp2/mesbah/tools/LiftOVer/grch38_to_grch37.over.chain R=/broad/hptmp/mesbah/RefSeq/human_g1k_v37.fasta VERBOSITY=DEBUG DISABLE_SORT=true

java -Xmx20G -jar ${Picard} LiftoverVcf \
	I=${sorted_VCF}.gz \
	O=${liftedVCF} \
	CHAIN=${hg38tohg37} \
	REJECT=${rejectVCF} \
	R=${refTarget} \
	VERBOSITY=DEBUG \
	DISABLE_SORT=true

# Need sorting to index!!
# tabix -f -p vcf ${outVCF}

######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

