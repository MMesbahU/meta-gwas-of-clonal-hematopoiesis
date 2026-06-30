#!/bin/bash

####################################################################
	## hg38
	# run time ~8hrs
	# qsub -wd /broad/hptmp/mesbah/RefSeq/dbSNP/tmpdir -R y -N prepare_dbSNP -l h_rt=20:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/11.etcs/add_rsids_dbSNP/run.dbSNP.v1.sh /broad/hptmp/mesbah/RefSeq/dbSNP/GCF_000001405.39.gz /broad/hptmp/mesbah/RefSeq/dbSNP/hg38 hg38 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/11.etcs/add_rsids_dbSNP/hg38.rename_chromosome.tsv
	
	## qsub -wd /broad/hptmp/mesbah/dbSNP/tmpdir -R y -N prepRSids -l h_rt=20:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/dbSNP/run.dbSNP.v1.sh /broad/hptmp/mesbah/dbSNP/GCF_000001405.39.gz /broad/hptmp/mesbah/dbSNP/hg38 hg38 /broad/hptmp/mesbah/dbSNP/hg38.rename_chromosome.tsv

	## hg37
	## qsub -wd /broad/hptmp/mesbah/dbSNP/tmpdir -R y -N prepRSids.hg37 -l h_rt=20:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/dbSNP/run.dbSNP.v1.sh /broad/hptmp/mesbah/dbSNP/GCF_000001405.25.gz /broad/hptmp/mesbah/dbSNP/hg37 hg37 /broad/hptmp/mesbah/dbSNP/hg37.rename_chromosome.tsv
	## for chr in {1..22}; do tabix GCF_000001405.25.gz -l | awk -v vars=${chr} 'NR==vars{print $1"\t"vars}' >> hg37.rename_chromosome.tsv; done; 
# echo -e "NC_000023.10\tX" >> hg37.rename_chromosome.tsv
#####################################################################
source /broad/software/scripts/useuse

use Tabix

# use Bcftools

######################################################################
################## Latest: May 25, 2021 (https://ftp.ncbi.nih.gov/snp/latest_release/VCF/: 2021-05-25)
###
## Note :
## BCFTOOLS bcftools 1.15.1 (Using htslib 1.15.1) does not work when contig info missing
## bcftools 1.10.2-57-g356d395 (Using htslib 1.10.2-82-g9de45b7)
###############
## hg38
# wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.39.{gz,gz.tbi}
# for chr in {1..23}; do tabix -h /broad/hptmp/mesbah/RefSeq/dbSNP/GCF_000001405.39.gz $(tabix /broad/hptmp/mesbah/RefSeq/dbSNP/GCF_000001405.39.gz -l | awk -v VARs=${chr} 'NR==VARs') | /medpop/esp2/mesbah/tools/bcftools/bcftools norm --multiallelics '-both' '-' | grep -v '^##' | awk -v CHR=${chr} '{print "chr"CHR":"$2":"$4":"$5"\tchr"CHR"\t"$2"\t"$3"\t"$4"\t"$5}' > /broad/hptmp/mesbah/RefSeq/dbSNP/hg38/chr${chr}.GCF_000001405.39.tsv && bgzip /broad/hptmp/mesbah/RefSeq/dbSNP/hg38/chr${chr}.GCF_000001405.39.tsv; done &
##
	### hg37
# wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.{gz,gz.tbi}

# for Chr in {1..23}; do tabix -h /broad/hptmp/mesbah/dbSNP/GCF_000001405.25.gz $(tabix /broad/hptmp/mesbah/dbSNP/GCF_000001405.25.gz -l | awk -v VARs=${Chr} 'NR==VARs') | bcftools norm --multiallelics '-both' '-' | grep -v '^##' | awk -v CHR=${Chr} '{print CHR":"$2":"$4":"$5"\tCHR\t"$2"\t"$3"\t"$4"\t"$5}' > /broad/hptmp/mesbah/dbSNP/hg37/chr${Chr}.GCF_000001405.37.tsv && bgzip /broad/hptmp/mesbah/dbSNP/hg37/chr${Chr}.GCF_000001405.39.tsv; done &

######################################################################

### Command line input
dbsnp_file=${1}

outDir=${2}

hgVersion=${3}

chr_rename=${4}

BCFTOOLs=/medpop/esp2/mesbah/tools/bcftools/bcftools
#################

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

