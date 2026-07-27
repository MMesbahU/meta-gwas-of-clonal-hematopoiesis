#!/bin/bash

## July 2, 2026
# AAF: [5%,1%,0.1%,0.01%,0.001%]
# qsub -t 1-22 -wd /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/rvas/tmpdir -R y -l h_vmem=20G -l h_rt=10:00:00 -pe smp 2 -binding linear:2 -N step2.mgbb_rvas.nochr /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/02.GWAS/3.RVAS/mgbb/3.step2.run_rvas.sh /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/rvas/step2 2 "0.05,0.01,0.001,0.0001,0.00001" /broad/hptmp/mesbah/dataset/mesbah_medpop_esp/GWAS_CHIP/MGBB53k/wes/pgen/MGB_53K_WES_genotype_variant_sample_QCed /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.Feb2024.tsv.gz "hasCHvaf02,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" "PC{1:10},Age_Genotyping,sqrAge_Genotyping" "Sex,ancestry_pred,Batch_CHIP_call" /broad/hptmp/mesbah/dataset/mesbah_medpop_esp/GWAS_CHIP/MGBB53k/MultiANC/step1/NULL_MODEL_pred.modified.list "/medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/rvas/mask_files/chr" ".MGB_53K_WES.vep.hc_LOF.annotation" "/medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/rvas/mask_files/nochr" ".MGB_53K_WES.vep.hc_LOF.setlist" "/medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/rvas/mask_files/lof.mask" "MGB_WES_RVAS" "MGB_53K_WES.RVAS_out_firth_aaf00001.nochr"

## Jun 30, 2026
# qsub -t 1-22 -wd /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/rvas/tmpdir -R y -l h_vmem=20G -l h_rt=10:00:00 -pe smp 2 -binding linear:2 -N step2.mgbb_rvas.nochr /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/02.GWAS/3.RVAS/mgbb/3.step2.run_rvas.sh /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/rvas/step2 2 "0.05,0.01,0.001,0.0001" /broad/hptmp/mesbah/dataset/mesbah_medpop_esp/GWAS_CHIP/MGBB53k/wes/pgen/MGB_53K_WES_genotype_variant_sample_QCed /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.Feb2024.tsv.gz "hasCHvaf02,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" "PC{1:10},Age_Genotyping,sqrAge_Genotyping" "Sex,ancestry_pred,Batch_CHIP_call" /broad/hptmp/mesbah/dataset/mesbah_medpop_esp/GWAS_CHIP/MGBB53k/MultiANC/step1/NULL_MODEL_pred.modified.list "/medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/rvas/mask_files/chr" ".MGB_53K_WES.vep.hc_LOF.annotation" "/medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/rvas/mask_files/nochr" ".MGB_53K_WES.vep.hc_LOF.setlist" "/medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/rvas/mask_files/lof.mask" "MGB_WES_RVAS" "MGB_53K_WES.RVAS_out_firth_aaf0001.nochr"

# Jun 29, 2026
# qsub -t 1-22 -wd /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/rvas/tmpdir -R y -l h_vmem=20G -l h_rt=10:00:00 -pe smp 2 -binding linear:2 -N step2.mgbb_rvas /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/02_2.RVAS/step2.run_rvas.sh /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/rvas/step2 2 "0.05,0.01,0.001,0.0001" /broad/hptmp/mesbah/dataset/mesbah_medpop_esp/GWAS_CHIP/MGBB53k/wes/pgen/MGB_53K_WES_genotype_variant_sample_QCed /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.Feb2024.tsv.gz "hasCHvaf02,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" "PC{1:10},Age_Genotyping,sqrAge_Genotyping" "Sex,ancestry_pred,Batch_CHIP_call" /broad/hptmp/mesbah/dataset/mesbah_medpop_esp/GWAS_CHIP/MGBB53k/MultiANC/step1/NULL_MODEL_pred.modified.list "/medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/02_2.RVAS/mask_files/chr" ".MGB_53K_WES.vep.hc_LOF.annotation" "/medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/02_2.RVAS/mask_files/chr" ".MGB_53K_WES.vep.hc_LOF.setlist" "/medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/02_2.RVAS/mask_files/lof.mask" "MGB_WES_RVAS" "MGB_53K_WES.RVAS_out_firth_aaf0001.chr"

###
# qsub -t 1-22 -wd /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/rvas/tmpdir -R y -l h_vmem=20G -l h_rt=10:00:00 -pe smp 2 -binding linear:2 -N step2.mgbb_rvas /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/02_2.RVAS/step2.run_rvas.sh /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/rvas/step2 2  

# for chr 1 2 {4..22}; do regenie --step 2 --pgen /broad/hptmp/mesbah/dataset/mesbah_medpop_esp/GWAS_CHIP/MGBB53k/wes/pgen/MGB_53K_WES_genotype_variant_sample_QCed_chr${chr} --covarFile /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.Feb2024.tsv.gz --phenoFile /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.Feb2024.tsv.gz --bt --firth --approx --phenoColList "hasCHvaf02,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" --covarColList "PC{1:10},Age_Genotyping,sqrAge_Genotyping" --catCovarList "Sex,ancestry_pred,Batch_CHIP_call" --pred /broad/hptmp/mesbah/dataset/mesbah_medpop_esp/GWAS_CHIP/MGBB53k/MultiANC/step1/NULL_MODEL_pred.modified.list --anno-file /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/02_2.RVAS/mask_files/chr${chr}.MGB_53K_WES.vep.hc_LOF.annotation --set-list /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/02_2.RVAS/mask_files/chr${chr}.MGB_53K_WES.vep.hc_LOF.setlist --mask-def /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/02_2.RVAS/mask_files/lof.mask --aaf-bins 0.05,0.01,0.001 --write-mask --bsize 200 --htp "MGB_WES_RVAS" --out /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/rvas/chr${chr}.MGB_53K_WES.RVAS_out_firth --vc-tests skato,acato-full --rgc-gene-p ; done &

# Jun 28, 2026

# regenie --step 2 --pgen /broad/hptmp/mesbah/dataset/mesbah_medpop_esp/GWAS_CHIP/MGBB53k/wes/pgen/MGB_53K_WES_genotype_variant_sample_QCed_chr3 --covarFile /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.Feb2024.tsv.gz --phenoFile /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.Feb2024.tsv.gz --bt --firth --approx --phenoColList "hasCHvaf02,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR" --covarColList "PC{1:10},Age_Genotyping,sqrAge_Genotyping" --catCovarList "Sex,ancestry_pred,Batch_CHIP_call" --pred /broad/hptmp/mesbah/dataset/mesbah_medpop_esp/GWAS_CHIP/MGBB53k/MultiANC/step1/NULL_MODEL_pred.list --anno-file /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/02_2.RVAS/mask_files/chr3.MGB_53K_WES.vep.hc_LOF.annotation --set-list /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/02_2.RVAS/mask_files/chr3.MGB_53K_WES.vep.hc_LOF.setlist --mask-def /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/02_2.RVAS/mask_files/lof.mask --aaf-bins 0.05,0.01,0.001 --write-mask --bsize 200 --out /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/chr3.MGB_53K_WES.RVAS_out_firth


##########################################################################
source /broad/software/scripts/useuse

source /home/unix/muddin/.bashrc

use Anaconda3

## to install
# conda create -n regenie_env -c conda-forge -c bioconda regenie
## Update from v3.1g to v3.1.3g to 4.1
# conda update -n regenie_env -c conda-forge -c bioconda regenie
# conda update -n base conda
# source activate regenie_env
# source activate /home/unix/muddin/micromamba/envs/regenie_env
# conda init bash

conda activate regenie_env

# conda deactivate
#########################################################################
chr=${SGE_TASK_ID}
outdir=${1}
cpus=${2} # 2
aaf_bin=${3} # 0.05,0.01,0.001,0.0001
pgen_prefix=${4} # /broad/hptmp/mesbah/dataset/mesbah_medpop_esp/GWAS_CHIP/MGBB53k/wes/pgen/MGB_53K_WES_genotype_variant_sample_QCed
pgen=${pgen_prefix}_chr${chr}
phenoFile=${5} # /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.Feb2024.tsv.gz
phenoColList=${6} # "hasCHvaf02,hasCHvaf10,hasDNMT3A,hasTET2,hasASXL1,hasSF,hasDDR"
covarColList=${7} # "PC{1:10},Age_Genotyping,sqrAge_Genotyping"
catCovarList=${8} # "Sex,ancestry_pred,Batch_CHIP_call"
pred=${9} # /broad/hptmp/mesbah/dataset/mesbah_medpop_esp/GWAS_CHIP/MGBB53k/MultiANC/step1/NULL_MODEL_pred.modified.list
annot_prefix=${10} # "/medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/02_2.RVAS/mask_files/chr"
annot_suffix=${11} # ".MGB_53K_WES.vep.hc_LOF.annotation"
anno_file=${annot_prefix}${chr}${annot_suffix}
set_prefix=${12} # "/medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/02_2.RVAS/mask_files/chr"
set_suffix=${13} # ".MGB_53K_WES.vep.hc_LOF.setlist"
setlist=${set_prefix}${chr}${set_suffix}
mask_def=${14} # "/medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/02_2.RVAS/mask_files/lof.mask"
cohort=${15} #"MGB_WES_RVAS"
outprefix=${16} # MGB_53K_WES.RVAS_out_firth.chr
outFile=${outdir}/${outprefix}${chr}
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



