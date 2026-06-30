#!/bin/bash

## for chr in {1..22}; do awk '$8>=0.3 && $6>=0.001 && $6<=0.999' /broad/ukbb/imputed_v3/ukb_mfi_chr${chr}_v3.txt > /broad/hptmp/mesbah/RefSeq/ukb200k/ukb_mfi_chr${chr}_v3.maf001_info30.txt; done &

# Samples to keep: zcat /medpop/esp/mesbah/GWAS_CHIP/inputUKB/aug11_2021_hg38/Chr1_additive_hasCHIP.regenie.ids.gz | awk 'NR>1' > /broad/hptmp/mesbah/ukb_chip/ukb200k_reference/ukb_200k_CHIP_GWAS.samples.tsv

# mkdir -p /medpop/esp2/mesbah/datasets/ukbb/ukbb_v3/tmpdir
# qsub -wd /medpop/esp2/mesbah/datasets/ukbb/ukbb_v3/tmpdir -R y -t 1-22 -l h_vmem=10G -l h_rt=20:00:00 -pe smp 6 -binding linear:6 -N prep_ukbb_bed /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.SNP_heritability/QuickPRS/prepUKB_imp_v3_GT.sh 6 /medpop/esp2/mesbah/datasets/ukbb/ukbb_v3/imp_v3
#
nNodes=${1}

# samples_2_keep=${2}
# snps_2_keep=${2}

outDir=${2}

chr=${SGE_TASK_ID}

snps_2_keep=/medpop/esp2/mesbah/datasets/ukbb/ukbb_v3/var/hapmap3.ukb_mfi_chr${chr}_v3.txt
###
######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
#########
/medpop/esp2/mesbah/tools/qcTool/bin/qctool_v2.0.7 \
	-g /broad/ukbb/imputed_v3/ukb_imp_chr${chr}_v3.bgen \
	-s /medpop/esp2/pradeep/UKBiobank/v3data/ukb7089_imp_chr3_v3_s487395.sample \
	-threads ${nNodes} \
	-incl-snpids ${snps_2_keep} \
	-og ${outDir}/ukb_imp_chr${chr}_v3 \
	-threshold 0.9 \
	-ofiletype binary_ped


######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

