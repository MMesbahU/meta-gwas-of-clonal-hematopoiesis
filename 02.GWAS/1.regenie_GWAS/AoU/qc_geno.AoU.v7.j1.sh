#!/bin/bash

# plink_files=gs://fc-aou-datasets-controlled/v7/microarray/plink/arrays.{bed,bim,fam}

# plink_docker=gcr.io/bick-aps2/plink2:polygenic

plink_path="$(dirname "${PLINK_FILES}")"
plink_base="$(basename "${PLINK_FILES//\*}")"
PLINK_Prefix=${plink_path}/${plink_base}
##
echo -e "all plink files ${PLINK_FILES}\n" 

echo -e "plink file prefix: ${PLINK_Prefix}\n"

### Plink filter
plink2 \
	--bfile ${PLINK_Prefix} \
	--maf 0.01 \
	--mac 100 \
	--geno 0.1 \
	--hwe 1e-50 \
	--mind 0.1 \
	--extract <(awk '$1~/chr[0-9]+/' ${PLINK_Prefix}.bim) \
	--write-snplist \
	--write-samples \
	--no-id-header \
	--threads $(nproc --all) \
	--out ${OUTPUT_PATH}/aou_qc_pass
##
# gsutil -m cp ${OUTPUT_PATH}/aou_qc_pass* gs/fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/aou_v7/qc_geno/

##
echo -e "qc is done \n"
# echo $(wc -l ${OUTPUT_PATH}/aou_qc_pass.snplist)

###
## Subset Genotype file
## keep only autosomes
# plink2 \
#	--bfile ${PLINK_Prefix} \
#	--extract ${OUTPUT_PATH}/aou_qc_pass.snplist \
#	--make-bed \
#	--out ${OUTPUT_PATH}/qc_arrays \
#	--threads $(nproc --all)

## 
# gsutil -m cp ${OUTPUT_PATH}/qc_arrays* gs/fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/aou_v7/qc_geno/

### Get top 10 PCs
## PCA Analysis
# plink2 \
#	--bfile ${OUTPUT_PATH}/qc_arrays \
#	--pca 10 \
#	--threads $(nproc --all) \
#	--out ${OUTPUT_PATH}/top10_PC_qc_arrays \
#	--memory 120000
	
## use e2-highmem-16 == 16 cpu 128 GB memory
###
# gsutil -m cp ${OUTPUT_PATH}/top10_PC_qc_arrays* gs/fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/aou_v7/qc_geno/


