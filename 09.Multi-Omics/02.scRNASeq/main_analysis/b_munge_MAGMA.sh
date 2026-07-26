#!/bin/bash

module load conda
conda activate scdrs

MAGMA_file='/data/cteu/2_users/zyu_lab/linke/chip_gwas_rev/Data/MAGMA_files/magma_chip_dnmt3a_tet2_merged.tsv'
out_file='/data/cteu/2_users/zyu_lab/linke/chip_gwas_rev/Data/scDRS/munged/chip_dnmt3a_tet2_munged_n1000.gs'
out_file_sensi="/data/cteu/2_users/zyu_lab/linke/chip_gwas_rev/Data/scDRS/munged/chip_dnmt3a_tet2_munged_n500.gs"

# Select top 1,000 genes and use z-score weights
scdrs munge-gs \
    --out-file ${out_file} \
    --zscore-file ${MAGMA_file} \
    --weight zscore \
    --n-max 1000


# for sensitivity
scdrs munge-gs \
    --out-file ${out_file_sensi} \
    --zscore-file ${MAGMA_file} \
    --weight zscore \
    --n-max 500
    
    
   
   
   
   
### ENSG files for ONEK1K 
MAGMA_file='/data/cteu/2_users/zyu_lab/linke/chip_gwas_rev/Data/MAGMA_files/magma_chip_dnmt3a_tet2_merged_ensg.tsv'
out_file='/data/cteu/2_users/zyu_lab/linke/chip_gwas_rev/Data/scDRS/munged/chip_dnmt3a_tet2_munged_n1000_ensg.gs'
out_file_sensi="/data/cteu/2_users/zyu_lab/linke/chip_gwas_rev/Data/scDRS/munged/chip_dnmt3a_tet2_munged_n500_ensg.gs"

# Select top 1,000 genes and use z-score weights
scdrs munge-gs \
    --out-file ${out_file} \
    --zscore-file ${MAGMA_file} \
    --weight zscore \
    --n-max 1000


# for sensitivity
scdrs munge-gs \
    --out-file ${out_file_sensi} \
    --zscore-file ${MAGMA_file} \
    --weight zscore \
    --n-max 500