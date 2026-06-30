#!/bin/sh
#$ -wd /medpop/esp2/btruong/Projects/logjobs
#$ -j y
#$ -pe smp 4 -R y -binding linear:4
#$ -l h_vmem=8G
#$ -l h_rt=20:00:00
#$ -N pops


source /broad/software/scripts/useuse
source ~/.my.bashrc
use GCC-5.2

reuse Python-3.6
reuse Anaconda
reuse Anaconda3


cd ${wdir}



python /medpop/esp2/btruong/Tools/pops/pops.py \
 --gene_annot_path /broad/hptmp/btruong/HDP/pops/gene_annot_jun10.txt \
 --feature_mat_prefix /broad/hptmp/btruong/HDP/pops/features_munged_1/pops_features \
 --num_feature_chunks 12 \
 --magma_prefix magma_${trait}_1 \
 --control_features_path /medpop/esp2/projects/software/PoPS/data/control.features \
 --verbose \
 --out_prefix pops_step2_${trait}
 
