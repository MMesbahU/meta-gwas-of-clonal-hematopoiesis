#!/bin/bash

## zcat /medpop/esp/mesbah/GWAS_CHIP/MGBB53k/chr1_22.MGBB53k.hasCHIP.tsv.gz | awk '(NR==1){print "SNPID"}(NR>1){print $1}' > /broad/hptmp/mesbah/RefSeq/mgbb53k/snpids.chr1_22.mgbb53k_ch_gwas.txt &

# Samples to keep: /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/multi_ancestry_summary/mgbb53k.multiAncestry.CHIP.ids

## qsub -wd /broad/hptmp/mesbah/RefSeq/tmpdir -R y -t 1-22 -l h_vmem=10G -l h_rt=20:00:00 -pe smp 6 -binding linear:6 -N prep_mgbb53k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis.v1/prepUKB/prepMGBB53k.Reference.sh 6 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/multi_ancestry_summary/mgbb53k.multiAncestry.CHIP.ids /broad/hptmp/mesbah/RefSeq/mgbb53k /broad/hptmp/mesbah/RefSeq/mgbb53k/snpids.chr1_22.mgbb53k_ch_gwas.txt

# mkdir -p /broad/hptmp/mesbah/RefSeq/mgbb53k
nNodes=${1}
samples_2_keep=${2}
outDir=${3}
chr=${SGE_TASK_ID}
snps_2_keep=${4}
###
######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
#########
# /medpop/esp2/mesbah/tools/qcTool/bin/qctool_v2.0.7 -g /medpop/esp2/projects/MGB_Biobank/imputation/53K_GSA/release/bgen/GSA_53K.merged.chr22.bgen -incl-samples /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/multi_ancestry_summary/mgbb53k.multiAncestry.CHIP.ids -incl-rsids snpids.chr1_22.mgbb53k_ch_gwas.txt -og test.chr22 -ofiletype binary_ped
#########
## Error with Sample file '/medpop/esp2/projects/MGB_Biobank/imputation/53K_GSA/release/bgen/GSA_53K.merged.chr${chr}.sample'
## works without it
/medpop/esp2/mesbah/tools/qcTool/bin/qctool_v2.0.7 \
	-g /medpop/esp2/projects/MGB_Biobank/imputation/53K_GSA/release/bgen/GSA_53K.merged.chr${chr}.bgen \
	-threads ${nNodes} \
	-incl-rsids ${snps_2_keep} \
	-incl-samples ${samples_2_keep} \
	-og ${outDir}/GSA_53K.merged.chr${chr}.mgbb53k_CHIP_samples \
	-ofiletype binary_ped

######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

