#!/bin/bash


### Tasks:
# 1. Prepare VCF
# 2. liftover hg38 to hg37
##

## qsub -R y -wd /broad/hptmp/mesbah/dataset/ch_gwas/tmpdir -t 1-$(wc -l /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/topmed_gwas_summary_list.list |awk '{print $1}') -pe smp 1 -binding linear:1 -l h_vmem=30G -l h_rt=10:00:00 -N hg38_to_hg37_topmed /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/11.etcs/topmed_hg38_to_hg37.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/topmed_gwas_summary_list.list /medpop/esp2/mesbah/tools/LiftOVer/grch38_to_grch37.over.chain /medpop/esp2/mesbah/tools/LiftOVer/picard.jar /broad/hptmp/mesbah/RefSeq/human_g1k_v37.fasta /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/vcf /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/vcf ".tsv.gz"

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
## 
# ls /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/*.regenie.tsv.gz| awk '{print $NF}' > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/topmed_gwas_summary_list.list
# zcat /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.EUR.TOPMed74k.chr1_22.hasCH.regenie.tsv.gz | head -2
# SNPID	RSID	REF	ALT	AAF	BETA	SE	P	N
# chr1:271618:C:T	rs11490246	T	C	0.000568786	-0.334440	0.494418	0.498767	40437

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
gwas_file_list=${1}

hg38tohg37=${2} # /medpop/esp2/mesbah/tools/LiftOVer/grch38_to_grch37.over.chain

Picard=${3} # /medpop/esp2/mesbah/tools/LiftOVer/picard.jar

refTarget=${4} #/broad/hptmp/mesbah/RefSeq/human_g1k_v37.fasta

hg38VCFDir=${5}

liftedOutDir=${6}

suffix_2_remove=${7} # ".tsv.gz"
###
## input files
gwas_file=$(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${gwas_file_list} )

inVCF=${hg38VCFDir}/hg38.$(basename $(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${gwas_file_list} ) ${suffix_2_remove}).vcf

sorted_VCF=${hg38VCFDir}/hg38.$(basename $(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${gwas_file_list} ) ${suffix_2_remove}).sorted.vcf

##
liftedVCF=${liftedOutDir}/hg37.$(basename $(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${gwas_file_list} ) ${suffix_2_remove}).vcf.gz

rejectVCF=${liftedOutDir}/rejected_hg37.$(basename $(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${gwas_file_list} ) ${suffix_2_remove}).vcf.gz

#############################
### Prepare VCF file
echo -e "input file: ${gwas_file} \n"

## convert to vcf
## https://genomics.ut.ee/en/tools/gwama-tutorial

# zcat /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.EUR.TOPMed74k.chr1_22.hasCH.regenie.tsv.gz | head -2
# SNPID RSID    REF     ALT     AAF     BETA    SE      P       N
# chr1:271618:C:T       rs11490246      T       C       0.000568786     -0.334440       0.494418        0.498767        40437
echo -e '##fileformat=VCFv4.2' > ${inVCF}

echo -e '##INFO=<ID=varID_hg38,Number=1,Type=String,Description="hg38 chr:pos:ref:alt ids">' >> ${inVCF}

echo -e '##INFO=<ID=RSID,Number=1,Type=String,Description="RSID">' >> ${inVCF}

echo -e '##INFO=<ID=EffectAllele,Number=1,Type=String,Description="Effect allele (Reference allele is the effect allele in topmed gwas)">' >> ${inVCF}

echo -e '##INFO=<ID=OtherAllele,Number=1,Type=String,Description="Alternate allele">' >> ${inVCF}

echo -e '##INFO=<ID=EAF,Number=A,Type=Float,Description="effect allele frequency (i.e. freq for REF allele). REF allele is the effect allele">' >> ${inVCF}

echo -e '##INFO=<ID=BETA,Number=A,Type=Float,Description="REGENIE effect estimate">' >> ${inVCF}

echo -e '##INFO=<ID=SE,Number=A,Type=Float,Description="REGENIE se of effect estimate">' >> ${inVCF}

echo -e '##INFO=<ID=Pvalue,Number=A,Type=String,Description="REGENIE P value">' >> ${inVCF}

echo -e '##INFO=<ID=N_Samples,Number=1,Type=Integer,Description="Number of samples with marker present">' >> ${inVCF}

echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' >> ${inVCF}

# zcat /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/topmed/GWAMA.EUR.TOPMed74k.chr1_22.hasCH.regenie.tsv.gz | head -2
# SNPID RSID    REF     ALT     AAF     BETA    SE      P       N
# chr1:271618:C:T       rs11490246      T       C       0.000568786     -0.334440       0.494418        0.498767        40437

zcat ${gwas_file} | awk 'BEGIN {OFS = "\t"} NR > 1 && $1 ~ /^chr[0-9]+/ {split($1, fields, ":"); chr = fields[1]; pos = fields[2]; ref = fields[3]; alt = fields[4]; varID = $1; rsid = $2; print chr, pos, rsid, ref, alt, ".", ".", "RSID=" rsid ";varID_hg38=" varID ";EffectAllele=" $4 ";OtherAllele=" $3 ";EAF=" $5 ";BETA=" $6 ";SE=" $7 ";Pvalue=" $8 ";N_Samples=" $9 }' >> ${inVCF}

## 
# Sort by chromosome
cat ${inVCF} | vcf-sort -c > ${sorted_VCF}

rm ${inVCF}

# compress  
gzip -f ${sorted_VCF}


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

rm ${sorted_VCF}.gz

######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

