#!/bin/bash

## hg37: qsub -t 1-22 -wd /broad/hptmp/mesbah/dbSNP/tmpdir -R y -N prepRSids.hg37 -l h_rt=20:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 /broad/hptmp/mesbah/dbSNP/run.dbSNP_hg37.sh /broad/hptmp/mesbah/dbSNP/GCF_000001405.25.gz /broad/hptmp/mesbah/dbSNP/hg37 hg37

source /broad/software/scripts/useuse

use Tabix

use Bcftools

# GRCh38
# wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.38.{gz,gz.tbi}
# for Chr in {1..22}; do tabix GCF_000001405.38.gz $(tabix GCF_000001405.38.gz -l | awk -v VARs=${Chr} 'NR==VARs') | awk -v CHR=${Chr} '{print "chr"CHR"\t"$2"\t"$3"\t"$4"\t"$5}' > grch38/Chr${Chr}.GCF_000001405.38.tsv && bgzip grch38/Chr${Chr}.GCF_000001405.38.tsv; done &

# GRCh37
# wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.{gz,gz.tbi}
# for chr in {1..22}; do tabix GCF_000001405.25.gz $(tabix GCF_000001405.25.gz -l | awk -v VARs=${chr} 'NR==VARs') | awk -v CHR=${chr} '{print "chr"CHR"\t"$2"\t"$3"\t"$4"\t"$5}' > grch37/chr${chr}.GCF_000001405.25.tsv && bgzip grch37/chr${chr}.GCF_000001405.25.tsv; done &

################## Latest: Aug 23, 2021
### Updated dbSNP

## hg38
# wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.39.{gz,gz.tbi}
# for chr in {1..23}; do tabix -h /broad/hptmp/mesbah/dbSNP/GCF_000001405.39.gz $(tabix /broad/hptmp/mesbah/dbSNP/GCF_000001405.39.gz -l | awk -v VARs=${chr} 'NR==VARs') | bcftools norm --multiallelics '-both' '-' | grep -v '^##' | awk -v CHR=${chr} '{print "chr"CHR":"$2":"$4":"$5"\tchr"CHR"\t"$2"\t"$3"\t"$4"\t"$5}' > /broad/hptmp/mesbah/dbSNP/hg38/chr${chr}.GCF_000001405.39.tsv && bgzip /broad/hptmp/mesbah/dbSNP/hg38/chr${chr}.GCF_000001405.39.tsv; done &
##
	### hg37
# wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.{gz,gz.tbi}

# for Chr in {1..23}; do tabix -h /broad/hptmp/mesbah/dbSNP/GCF_000001405.25.gz $(tabix /broad/hptmp/mesbah/dbSNP/GCF_000001405.25.gz -l | awk -v VARs=${Chr} 'NR==VARs') | bcftools norm --multiallelics '-both' '-' | grep -v '^##' | awk -v CHR=${Chr} '{print CHR":"$2":"$4":"$5"\tCHR\t"$2"\t"$3"\t"$4"\t"$5}' > /broad/hptmp/mesbah/dbSNP/hg37/chr${Chr}.GCF_000001405.37.tsv && bgzip /broad/hptmp/mesbah/dbSNP/hg37/chr${Chr}.GCF_000001405.39.tsv; done &

### 
dbsnp_file=${1}

outDir=${2}

hgVersion=${3}

chr=${SGE_TASK_ID}

## HG37
tabix -h ${dbsnp_file} $(tabix ${dbsnp_file} -l | awk -v VARs=${chr} 'NR==VARs') | bcftools norm --multiallelics '-both' '-' | grep -v '^##' | awk -v CHR=${chr} '{print CHR":"$2":"$4":"$5"\t"CHR"\t"$2"\t"$3"\t"$4"\t"$5}' | bgzip -c > ${outDir}/chr${chr}.dbSNP.${hgVersion}.tsv.gz

