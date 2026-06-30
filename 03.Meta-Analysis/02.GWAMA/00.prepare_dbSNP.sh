#!/bin/bash

####################################################################
	## hg38
	# run time ~8hrs
	# qsub -wd /broad/hptmp/mesbah/RefSeq/dbSNP/tmpdir -R y -N prepare_dbSNP -l h_rt=20:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/00.prepare_dbSNP.sh /medpop/esp2/mesbah/tools/bcftools/bcftools hg38 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/11.etcs/add_rsids_dbSNP/hg38.rename_chromosome.tsv
##
# Task performed:
# 1. download dbSNP hg38
# 2. Normalize multiallelic SNP and INDELs (bcftools 'norm' function); prepare index file
# 3. Change contig names to chr1-chrX; prepare index file
### Run time ~8 hrs
#####################################################################
source /broad/software/scripts/useuse

use Tabix

# use Bcftools

### Command line input
BCFTOOLs=${1} # /medpop/esp2/mesbah/tools/bcftools/bcftools

hgVersion=${2} # hg38

chr_rename=${3}
##
mkdir -p /broad/hptmp/mesbah/RefSeq/dbSNP/hg38

outDir=/broad/hptmp/mesbah/RefSeq/dbSNP/hg38

### Download dbSNP hg38
wget -P /broad/hptmp/mesbah/RefSeq/dbSNP https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.39.{gz,gz.tbi}

dbsnp_file=/broad/hptmp/mesbah/RefSeq/dbSNP/GCF_000001405.39.gz

######################################################################
################## Latest: May 25, 2021 (https://ftp.ncbi.nih.gov/snp/latest_release/VCF/: 2021-05-25)
###
## Note :
## BCFTOOLS bcftools 1.15.1 (Using htslib 1.15.1) does not work when contig info missing
## bcftools 1.10.2-57-g356d395 (Using htslib 1.10.2-82-g9de45b7)
###############
## hg38
# wget -P /broad/hptmp/mesbah/RefSeq/dbSNP https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.39.{gz,gz.tbi}
##
### hg37
# wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.{gz,gz.tbi}
######################################################################

######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
#########

#### Normalize multiallelic variants 
tabix \
	-h ${dbsnp_file} $(tabix ${dbsnp_file} -l | head -23) | \
	${BCFTOOLs} norm \
		--multiallelics '-both' '-' \
		-o ${outDir}/${hgVersion}.$(basename ${dbsnp_file} ".gz").vcf.gz \
		-Oz
## Prep index 
tabix -p vcf ${outDir}/${hgVersion}.$(basename ${dbsnp_file} ".gz").vcf.gz

## Change chromosome name
${BCFTOOLs} annotate ${outDir}/${hgVersion}.$(basename ${dbsnp_file} ".gz").vcf.gz \
	--rename-chrs ${chr_rename} \
	-Oz \
	-o ${outDir}/chr1_X.${hgVersion}.$(basename ${dbsnp_file} ".gz").vcf.gz

## Prep index
tabix -p vcf ${outDir}/chr1_X.${hgVersion}.$(basename ${dbsnp_file} ".gz").vcf.gz

## Remove intermediate file
rm ${outDir}/${hgVersion}.$(basename ${dbsnp_file} ".gz").vcf.*

# bcftools annotate -a /broad/hptmp/mesbah/dbSNP/hg38/chr1.dbSNP.hg38.vcf.gz -c ID /broad/hptmp/mesbah/ukb_chip/eurUKB/lifted.ukb_Chr1_22.hasCHIP.regenie.vcf.gz -o /broad/hptmp/mesbah/ukb_chip/eurUKB/dbSNP.lifted.ukb_Chr1_22.hasCHIP.regenie.vcf.gz -O z

######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

