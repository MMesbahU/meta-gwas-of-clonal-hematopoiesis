#!/bin/bash

source /broad/software/scripts/useuse

use .bgen-1.1.4

QCTool_v2="/medpop/esp2/mesbah/tools/qcTool/bin/qctool_v2.0.7"

LDstore2="/medpop/esp2/mesbah/tools/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64"

FINEMAP="/medpop/esp2/mesbah/tools/finemap_v1.4.2_x86_64/finemap_v1.4.2_x86_64"
# EAF 1%
# qsub -t 1-$(wc -l /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/finemap_loci_500kb.hasCH.MultiANC.multianc_cohort6.txt |awk '{print $1}') -wd /broad/hptmp/mesbah/RefSeq/tmpdir -R y -l h_vmem=20G -l h_rt=20:00:00 -pe smp 4 -binding linear:4 -N MGBB_CHIP_MultiANC.mgbb53k.bgen_ldstore_finemap_eaf01 /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/07.Fine-mapping/04.bgen_to_ldstore_to_finemap.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/finemap_loci_500kb.hasCH.MultiANC.multianc_cohort6.txt CHIP_MultiANC /medpop/esp2/projects/MGB_Biobank/imputation/53K_GSA/release/bgen /medpop/esp2/mesbah/projects/Meta_GWAS/n650k/multi_ancestry_summary/mgbb53k.multiAncestry.CHIP.ids /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/mgbb /broad/hptmp/mesbah/RefSeq/mgbb53k/region hasCH.chr1_22.MultiANC.0.01_0.99 50305 /broad/hptmp/mesbah/RefSeq/mgbb53k/ldstore 4 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/output

# EAF 0.1%
# qsub -t 1-$(wc -l /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/finemap_loci_500kb.hasCH.MultiANC.multianc_cohort6.txt |awk '{print $1}') -wd /broad/hptmp/mesbah/RefSeq/tmpdir -R y -l h_vmem=20G -l h_rt=20:00:00 -pe smp 4 -binding linear:4 -N CHIP_MultiANC.mgbb53k.bgen_ldstore_finemap_eaf001 /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/07.Fine-mapping/04.bgen_to_ldstore_to_finemap.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/finemap_loci_500kb.hasCH.MultiANC.multianc_cohort6.txt MGBB_CHIP_MultiANC /medpop/esp2/projects/MGB_Biobank/imputation/53K_GSA/release/bgen /medpop/esp2/mesbah/projects/Meta_GWAS/n650k/multi_ancestry_summary/mgbb53k.multiAncestry.CHIP.ids /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/mgbb /broad/hptmp/mesbah/RefSeq/mgbb53k/region hasCH.chr1_22.MultiANC.0.001_0.999 50305 /broad/hptmp/mesbah/RefSeq/mgbb53k/ldstore 4 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/output

##
LOCI=${1} # locus file

trait=${2}   # e.g. CHIP_multianc UKB, AoU, MGBB

BGEN_DIR=${3} # /medpop/esp2/projects/MGB_Biobank/imputation/53K_GSA/release/bgen

SAMPLE=${4} # /medpop/esp2/mesbah/projects/Meta_GWAS/n650k/multi_ancestry_summary/mgbb53k.multiAncestry.CHIP.ids

inputDir=${5} # /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/mgbb

bgenOUTDIR=${6} # 

gwas_suffix=${7} # hasCH.chr1_22.MultiANC

sampleN=${8} # 50305

ldstore_outdir=${9} #/broad/hptmp/mesbah/RefSeq/mgbb53k/ldstore

cpus=${10}

finemap_outdir_root=${11} # /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/output

mkdir -p ${finemap_outdir_root}/${trait}

finemap_outdir=${finemap_outdir_root}/${trait}

row_number=${SGE_TASK_ID}

######
while read CHR START END LEAD; do

    LOCUS=${LEAD}_${CHR}_${START}_${END}
    
    BGEN_file=${BGEN_DIR}/GSA_53K.merged.${CHR}.bgen

    RSID_file=${inputDir}/rsid.${CHR}_${START}_${END}_${LEAD}.${gwas_suffix}.z.txt
    
    gwas_Z=${inputDir}/${CHR}_${START}_${END}_${LEAD}.${gwas_suffix}.z
    
    out_bgen=${bgenOUTDIR}/${CHR}_${START}_${END}.${LEAD}.${gwas_suffix}.bgen
    
    qcd_bgen=${bgenOUTDIR}/$(basename ${out_bgen} ".bgen").MGBB_sample_N${sampleN}.bgen
    
    qcd_bgen_sample=${bgenOUTDIR}/$(basename ${out_bgen} ".bgen").MGBB_sample_N${sampleN}.sample


    echo -e "Running locus: ${LOCUS}"

    #############################################
    # STEP 1: Extract region from BGEN
    #############################################
# source /broad/software/scripts/useuse

# use .bgen-1.1.4
#
bgenix \
	-g ${BGEN_file} \
	-incl-rsids ${RSID_file} > ${out_bgen}
#
bgenix -g ${out_bgen} -index

## QCTool
## QCTool_v2="/medpop/esp2/mesbah/tools/qcTool/bin/qctool_v2.0.7"
${QCTool_v2} \
	-g ${out_bgen} \
	-incl-samples ${SAMPLE} \
	-og ${qcd_bgen} \
	-os ${qcd_bgen_sample} \
	-ofiletype bgen_v1.2 \
	-threads ${cpus}
#
bgenix -g ${qcd_bgen} -index

    #############################################
    # STEP 2: LDSTORE2 (recommended)
    #############################################
# ldstore z file
# rsid chromosome position allele1 allele2
bgen=${qcd_bgen}

bgi=${qcd_bgen}.bgi

n_samples=${sampleN}

ldstore_Z=${ldstore_outdir}/ldstoreZ.$(basename ${gwas_Z})

ldstore_master=${ldstore_outdir}/ldstore.$(basename ${gwas_Z} ".z").master

ldstore_bdose=${ldstore_outdir}/ldstore.$(basename ${gwas_Z} ".z").bdose

ldstore_bcor=${ldstore_outdir}/ldstore.$(basename ${gwas_Z} ".z").bcor

ldstore_ld=${ldstore_outdir}/ldstore.$(basename ${gwas_Z} ".z").ld

## LDstore2="/medpop/esp2/mesbah/tools/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64" 
##

cut -d ' ' -f1-5 ${gwas_Z} > ${ldstore_Z}

## Master file LDstore
cat > ${ldstore_master} << EOF
z;bgen;bgi;bcor;bdose;n_samples
${ldstore_Z};${bgen};${bgi};${ldstore_bcor};${ldstore_bdose};${n_samples}
EOF

## --n-threads ${thr} --memory ${mem} 
# Error : '--write-bcor' and '--bcor-to-text' cannot be used together! --bcor-to-text
#
${LDstore2} \
	--in-files ${ldstore_master} \
        --write-bdose \
        --write-bcor \
	--n-threads ${cpus}
###

###
finemap_master=${finemap_outdir}/finemap.$(basename ${gwas_Z} ".z").master

finemap_snp=${finemap_outdir}/finemap.$(basename ${gwas_Z} ".z").snp

finemap_log=${finemap_outdir}/finemap.$(basename ${gwas_Z} ".z").log

finemap_config=${finemap_outdir}/finemap.$(basename ${gwas_Z} ".z").config

finemap_cred=${finemap_outdir}/finemap.$(basename ${gwas_Z} ".z").cred
    
#############################################
    # STEP 5: Create FINEMAP master file
    #############################################

    cat > ${finemap_master} <<EOF
z;bcor;snp;config;cred;log;n_samples
${gwas_Z};${ldstore_bcor};${finemap_snp};${finemap_config};${finemap_cred};${finemap_log};${n_samples}
EOF

    #############################################
    # STEP 6: Run FINEMAP
    #############################################

${FINEMAP} \
        --sss \
        --in-files ${finemap_master} \
        --n-threads ${cpus} \
        --n-causal-snps 10 \
        --corr-config 0.95 \
        --log \
        --dataset 1

##### Rename finemap files
mv ${gwas_Z} ${gwas_Z}.txt
mv ${finemap_snp} ${finemap_snp}.txt
mv ${finemap_config} ${finemap_config}.txt
# mv ${finemap_cred} ${finemap_cred}.txt
# mv ${finemap_log} ${finemap_log}.txt

# remove begn files
rm ${ldstore_Z} ${bgen} ${bgi} ${ldstore_bcor} ${ldstore_bdose}
#####

done < <(cut -f1-4 ${LOCI} | awk -v row_num=${row_number} 'NR==row_num{print $0}' )

