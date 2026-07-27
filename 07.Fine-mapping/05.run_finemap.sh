#!/bin/bash


FINEMAP="/medpop/esp2/mesbah/tools/finemap_v1.4.2_x86_64/finemap_v1.4.2_x86_64"

# qsub -t 1-$(wc -l /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/finemap_loci_500kb.hasCH.MultiANC.multianc_cohort6.txt |awk '{print $1}') -wd /broad/hptmp/mesbah/RefSeq/tmpdir -R y -l h_vmem=20G -l h_rt=20:00:00 -pe smp 4 -binding linear:4 -N finemap.chip.mgbb53k /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/07.Fine-mapping/05.run_finemap.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/finemap_loci_500kb.hasCH.MultiANC.multianc_cohort6.txt MGBB /medpop/esp2/projects/MGB_Biobank/imputation/53K_GSA/release/bgen /medpop/esp2/mesbah/projects/Meta_GWAS/n650k/multi_ancestry_summary/mgbb53k.multiAncestry.CHIP.ids /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/mgbb /broad/hptmp/mesbah/RefSeq/mgbb53k/region hasCH.chr1_22.MultiANC 50305 /broad/hptmp/mesbah/RefSeq/mgbb53k/ldstore 4 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/output


##
LOCI=${1} # locus file

COHORT=${2}   # e.g. UKB, AoU, MGBB

BGEN_DIR=${3} # /medpop/esp2/projects/MGB_Biobank/imputation/53K_GSA/release/bgen

SAMPLE=${4} # /medpop/esp2/mesbah/projects/Meta_GWAS/n650k/multi_ancestry_summary/mgbb53k.multiAncestry.CHIP.ids

inputDir=${5} # /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/mgbb

bgenOUTDIR=${6} # 

gwas_suffix=${7} # hasCH.chr1_22.MultiANC

sampleN=${8} # 50305

ldstore_outdir=${9} #/broad/hptmp/mesbah/RefSeq/mgbb53k/ldstore

cpus=${10}

finemap_outdir=${11} # /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/output

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

    bgen=${qcd_bgen}

    bgi=${qcd_bgen}.bgi

    n_samples=${sampleN}

    ldstore_Z=${ldstore_outdir}/ldstoreZ.$(basename ${gwas_Z})

    ldstore_master=${ldstore_outdir}/ldstore.$(basename ${gwas_Z} ".z").master

    ldstore_bdose=${ldstore_outdir}/ldstore.$(basename ${gwas_Z} ".z").bdose

    ldstore_bcor=${ldstore_outdir}/ldstore.$(basename ${gwas_Z} ".z").bcor

    ldstore_ld=${ldstore_outdir}/ldstore.$(basename ${gwas_Z} ".z").ld

    finemap_master=${finemap_outdir}/finemap.$(basename ${gwas_Z} ".z").master

    finemap_snp=${finemap_outdir}/finemap.$(basename ${gwas_Z} ".z").snp

    finemap_log=${finemap_outdir}/finemap.$(basename ${gwas_Z} ".z").log

    finemap_config=${finemap_outdir}/finemap.$(basename ${gwas_Z} ".z").config

    finemap_cred=${finemap_outdir}/finemap.$(basename ${gwas_Z} ".z").cred
    
    echo -e "Running locus: ${LOCUS}"

    #############################################
    # STEP 4: Match LD variants with GWAS variants
    #############################################

 #   awk '{print $1}' ${WORK}/ld.bcor.vars > ${WORK}/ld.snps

 #   head -1 ${WORK}/z_raw.txt > ${WORK}/z.tmp

 #   grep -wFf ${WORK}/ld.snps ${WORK}/z_raw.txt >> ${WORK}/z.tmp

 #   mv ${WORK}/z.tmp ${WORK}/region.z

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


done < <(cut -f1-4 ${LOCI} | awk -v row_num=${row_number} 'NR==row_num{print $0}' )


