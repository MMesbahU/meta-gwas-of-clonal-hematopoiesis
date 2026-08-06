#!/bin/bash

# for chr in {1..22} X; do awk -v chrom=${chr} '{print $1"\t"chrom"\t"$3"\t"$4}' ~/exome_v7/mask_files/chr${chr}.AoUv7_exome.genotypeQCed.variantQCed.sampleQCed.vep.hc_LOF.setlist > ~/exome_v7/mask_files/nochr${chr}.AoUv7_exome.genotypeQCed.variantQCed.sampleQCed.vep.hc_LOF.setlist; done 

# Jun 29, 2026
# chr=22; bash /home/jupyter/meta-gwas-of-clonal-hematopoiesis/02.GWAS/3.RVAS/aou/3.step2.run_rvas.sh $chr /home/jupyter/workspace/rw-migration-aou-rw-424d5361/aou_v6/aou_v7/rvas/step2 4 "0.05,0.01,0.001,0.0001" /home/jupyter/exome_v7/genotypeQCed/AoUv7_exome.genotypeQCed /home/jupyter/workspace/rw-migration-aou-rw-424d5361/aou_v6/my_CH/phenoCH_AoU250k.fid0_iid.qcd_myeloidCA_rel_NA.26Aug2023.tsv "hasCH,hasCHvaf10,hasDTA,hasDDR,hasSF,hasDNMT3A,hasTET2,hasASXL1,hasPPM1D,hasTP53,hasSF3B1" Age_biosample_collection,sqrAge_biosample_collection,PC{1:10} dragen_sex_ploidy,ancestry_pred_other,site_id,Batch_CH /home/jupyter/workspace/rw-migration-aou-rw-424d5361/aou_v6/aou_v7/gwas/AoU_multiAnc_2023/step1/NULL_MODEL_pred_modified.verily_vm.list /home/jupyter/workspace/rw-migration-aou-rw-424d5361/aou_v6/aou_v7/rvas/mask_files/chr .AoUv7_exome.genotypeQCed.variantQCed.sampleQCed.vep.hc_LOF.annotation /home/jupyter/workspace/rw-migration-aou-rw-424d5361/aou_v6/aou_v7/rvas/mask_files/nochr .AoUv7_exome.genotypeQCed.variantQCed.sampleQCed.vep.hc_LOF.setlist /home/jupyter/workspace/rw-migration-aou-rw-424d5361/aou_v6/aou_v7/rvas/mask_files/lof.mask AoU_WGS_RVAS AoU_v7.250k_WGS.RVAS_out_firth.chr
#
#

##########################################################################
#########################################################################
# chr=${SGE_TASK_ID}
chr=${1}
outdir=${2}
cpus=${3} # 2
aaf_bin=${4} # 0.05,0.01,0.001,0.0001
pgen_prefix=${5} # /broad/hptmp/mesbah/dataset/mesbah_medpop_esp/GWAS_CHIP/MGBB53k/wes/pgen/MGB_53K_WES_genotype_variant_sample_QCed
pgen=${pgen_prefix}.chr${chr}
phenoFile=${6} # /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.Feb2024.tsv.gz
phenoColList=${7} # "hasCHvaf02,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR"
covarColList=${8} # "PC{1:10},Age_Genotyping,sqrAge_Genotyping"
catCovarList=${9} # "Sex,ancestry_pred,Batch_CHIP_call"
pred=${10} # /broad/hptmp/mesbah/dataset/mesbah_medpop_esp/GWAS_CHIP/MGBB53k/MultiANC/step1/NULL_MODEL_pred.modified.list
annot_prefix=${11} # "/medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/02_2.RVAS/mask_files/chr"
annot_suffix=${12} # ".MGB_53K_WES.vep.hc_LOF.annotation"
anno_file=${annot_prefix}${chr}${annot_suffix}
set_prefix=${13} # "/medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/02_2.RVAS/mask_files/chr"
set_suffix=${14} # ".MGB_53K_WES.vep.hc_LOF.setlist"
setlist=${set_prefix}${chr}${set_suffix}
mask_def=${15} # "/medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/02_2.RVAS/mask_files/lof.mask"
cohort=${16} #"MGB_WES_RVAS"
outprefix=${17} # MGB_53K_WES.RVAS_out_firth.chr
outFile=${outdir}/${outprefix}${chr}

# pgen=${pgen_prefix}.chr${chr}
# outFile=${outdir}/${outprefix}${chr}
##
#########################################################################
######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
#########
echo -e "##########################\n"
echo -e "Job interval: chr${chr}\t\n"
echo -e "##########################\n"
########
## Run Regenie
regenie \
	--step 2 \
	--pgen ${pgen} \
	--covarFile ${phenoFile} \
	--phenoFile ${phenoFile} \
	--bt \
	--firth \
	--approx \
	--phenoColList ${phenoColList} \
	--covarColList ${covarColList} \
	--catCovarList ${catCovarList} \
	--pred ${pred} \
	--anno-file ${anno_file} \
	--set-list ${setlist} \
	--mask-def ${mask_def} \
	--aaf-bins ${aaf_bin} \
	--write-mask \
	--bsize 200 \
	--htp ${cohort} \
	--out ${outFile} \
	--vc-tests skato,acato-full \
	--rgc-gene-p \
	--threads ${cpus}


######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########



