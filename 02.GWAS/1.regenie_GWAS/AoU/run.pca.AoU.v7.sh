#!/bin/bash
#######################
# Author: Md Mesbah Uddin <mdmesbah@gmail.com>
# August 10, 2023
#######################

### Run: bash run.pca.AoU.v7.sh pca.AoU.v7.sh QC_2023Aug11 gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/aou_v7/qc_geno/QC_2023Aug11/aou_qc_pass.snplist

### cp scripts to: 
PROJECT_PATH="gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/aou_v7/qc_geno"

plink_files="gs://fc-aou-datasets-controlled/v7/microarray/plink/arrays*"

PLINK_DOCKER="gcr.io/bick-aps2/plink2:polygenic"


## use e2-highmem-16 == 16 cpu 128 GB memory
e2_highmem_16='e2-highmem-16' ## or high cpu: e2-highcpu-32

##
PLINK_SCRIPT=${1}

job_name=${2}

snp_list=${3} 
# Launch step 1 and block until completion
echo -e "Launching plinks for approximate PCA\n"
###
#local
AOU_NETWORK=network
# local
AOU_SUBNETWORK=subnetwork

DSUB_USER_NAME="$(echo "${OWNER_EMAIL}" | cut -d@ -f1)"

## Dsub 
# Enable exit on error
# set -o errexit
## https://cloud.google.com/compute/vm-instance-pricing
dsub \
        --user-project ${GOOGLE_PROJECT} \
        --project ${GOOGLE_PROJECT} \
        --network ${AOU_NETWORK} \
        --subnetwork ${AOU_SUBNETWORK} \
        --service-account $(gcloud config get-value account) \
        --user ${DSUB_USER_NAME} \
        --provider google-cls-v2 \
        --regions us-central1 \
        --disk-type pd-standard \
        --disk-size 200 \
        --machine-type=${e2_highmem_16} \
        --image ${PLINK_DOCKER} \
        --logging ${PROJECT_PATH}/${job_name}/dsub-logs \
        --input PLINK_FILES=${plink_files} \
	--input SNP_LIST=${snp_list} \
        --output-recursive OUTPUT_PATH=${PROJECT_PATH}/${job_name} \
        --timeout '1w' \
        --name ${job_name}-pca \
        --label 'step=0' \
        --label 'script=aprox_pca_qc_snp' \
        --script ${PLINK_SCRIPT}

