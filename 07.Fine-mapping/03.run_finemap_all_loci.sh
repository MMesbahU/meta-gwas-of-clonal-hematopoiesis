#!/bin/bash

## bash run_finemap_all_loci.sh finemap_loci.txt UKB /data/ukb/bgen /data/ukb/sample/ukb.sample /data/ukb/gwas/trait.gz finemap_output


set -euo pipefail

LOCI=$1
COHORT=$2   # e.g. UKB, AoU, MGBB

BGEN_DIR=$3
SAMPLE=$4
GWAS=$5
OUTDIR=$6

N=420000   # adjust per cohort

mkdir -p ${OUTDIR}

while read CHR START END LEAD; do

    LOCUS=${LEAD}_${CHR}_${START}_${END}
    WORK=${OUTDIR}/${LOCUS}

    mkdir -p ${WORK}

    echo "Running locus: ${LOCUS}"

    #############################################
    # STEP 1: Extract region from BGEN
    #############################################

    bgenix \
        -g ${BGEN_DIR}/chr${CHR}.bgen \
        -incl-range ${CHR}:${START}-${END} \
        > ${WORK}/region.bgen

    bgenix \
        -g ${WORK}/region.bgen \
        -index

    #############################################
    # STEP 2: LDSTORE2 (recommended)
    #############################################

    ldstore_v2 \
        --bgen ${WORK}/region.bgen \
        --sample ${SAMPLE} \
        --write-bcor \
        --out ${WORK}/ld

    #############################################
    # STEP 3: Build FINEMAP Z file
    #############################################

    zcat ${GWAS} | awk -v chr=$CHR -v start=$START -v end=$END '
    BEGIN{
        OFS=" "
        print "rsid","chromosome","position","allele1","allele2","maf","beta","se"
    }
    NR==1{
        for(i=1;i<=NF;i++){
            if($i=="rs_number") id=i
            if($i=="reference_allele") ref=i
            if($i=="other_allele") alt=i
            if($i=="eaf") maf=i
            if($i=="beta") b=i
            if($i=="se") se=i
        }
    }
    NR>1{
        split($id,a,":")
        c=a[1]
        p=a[2]

        if(c==chr && p>=start && p<=end)
            print $id,c,p,$alt,$ref,$maf,$b,$se
    }' > ${WORK}/z_raw.txt

    #############################################
    # STEP 4: Match LD variants with GWAS variants
    #############################################

    awk '{print $1}' ${WORK}/ld.bcor.vars > ${WORK}/ld.snps

    head -1 ${WORK}/z_raw.txt > ${WORK}/z.tmp

    grep -wFf ${WORK}/ld.snps ${WORK}/z_raw.txt >> ${WORK}/z.tmp

    mv ${WORK}/z.tmp ${WORK}/region.z

    #############################################
    # STEP 5: Create FINEMAP master file
    #############################################

    cat > ${WORK}/master <<EOF
z;bcor;snp;config;cred;log;n_samples
region.z;ld.bcor;region.snp;region.config;region.cred;region.log;${N}
EOF

    #############################################
    # STEP 6: Run FINEMAP
    #############################################

    finemap_v1.4_x86_64 \
        --sss \
        --in-files ${WORK}/master \
        --n-causal-snps 10 \
        --corr-config 0.95 \
        --log

done < ${LOCI}
