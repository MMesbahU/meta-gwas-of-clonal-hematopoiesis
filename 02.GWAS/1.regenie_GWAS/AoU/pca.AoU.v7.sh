#!/bin/bash

# plink_files=gs://fc-aou-datasets-controlled/v7/microarray/plink/arrays.{bed,bim,fam}

# plink_docker=gcr.io/bick-aps2/plink2:polygenic

plink_path="$(dirname "${PLINK_FILES}")"
plink_base="$(basename "${PLINK_FILES//\*}")"
PLINK_Prefix=${plink_path}/${plink_base}

### Get top 20 PCs
### use Approximate approach
## exclude indels for PCA
## PCA Analysis
plink2 \
	--bfile ${PLINK_Prefix} \
	--extract <( sed 's:\::\t:g' ${SNP_LIST} | awk '{snp=$1":"$2":"$3":"$4;REF=$3;ALT=$4}( (REF=="A"||REF=="C"||REF=="G"||REF=="T") && (ALT=="A"||ALT=="C"||ALT=="G"||ALT=="T")){print snp}' ) \
	--pca 20 approx \
	--threads $(nproc --all) \
	--out ${OUTPUT_PATH}/top20_approx_PC_qc_arrays_SNPs \
	--memory 120000
	
## use e2-highmem-16 == 16 cpu 128 GB memory
###


