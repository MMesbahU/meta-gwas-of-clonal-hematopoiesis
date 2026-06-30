#!/bin/bash

# while read ANC; do qsub -t 1-$(wc -l /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/prs/list.ldak_input.txt | awk '{print $1}') -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/tmpdir -N ${ANC}.chip_megaprs -l h_rt=4:00:00 -l h_vmem=20G -pe smp 4 -binding linear:4 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.SNP_heritability/QuickPRS/get_bayesR.sh /medpop/esp2/mesbah/tools/LDAK5.2/ldak5.2.linux /medpop/esp2/mesbah/tools/LDAK6.0/ldak6.linux 4 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/prs/list.ldak_input.txt /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/prs/per_snp/${ANC} /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/prs/mega/${ANC} /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/${ANC}.hapmap/${ANC}.hapmap.bld.ldak.quickprs.tagging /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/${ANC}.hapmap/${ANC}.hapmap.bld.ldak.quickprs.matrix /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/${ANC}.hapmap/${ANC}.hapmap /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/${ANC}.hapmap/highld.snps; done < <(echo -e "afr\ngbr")

LDAK5=${1} # /medpop/esp2/mesbah/tools/LDAK5.2/ldak5.2.linux

LDAK6=${2}

maxThreads=${3}

##
list_hg37_summary_file=${4} # /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/prs/list.ldak_input.txt

ldak_inpuDir_persnp=${5}

ldak_inpuDir_megaprs=${6}

tagfile=${7} # /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/gbr.hapmap/gbr.hapmap.bld.ldak.quickprs.tagging 

matrix=${8} # /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/gbr.hapmap/gbr.hapmap.bld.ldak.quickprs.matrix

ref_cors=${9} # /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/gbr.hapmap/gbr.hapmap

ref_highld=${10} # /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/gbr.hapmap/highld.snps
####
### get GWAS file from t array file list
gwas_summary_ldak_format=$(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${list_hg37_summary_file} )

per_snp_hers=${ldak_inpuDir_persnp}/hers_ldak.$(basename $(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${list_hg37_summary_file} ) ".txt")

mega_prs_weight=${ldak_inpuDir_megaprs}/mega_prs.$(basename $(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${list_hg37_summary_file} ) ".txt")
####

## Step 0:
##### Prepare LDAK formated GWAS summary
# Run: /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.SNP_heritability/QuickPRS/prep_LDAK_summary.sh 
# ls /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/prs/ldak.hg37.GWAMA.chr1_22.has*.txt| awk '{print $NF}' | sort -V > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/prs/list.ldak_input.txt

# qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/tmpdir -N prep_ldak_sum -l h_rt=1:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 -t 1-$(wc -l /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/prs/list_hg37_aou_top_mgbb.noukb.txt | awk '{print $1}') /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.SNP_heritability/QuickPRS/prep_LDAK_summary.sh /medpop/esp2/mesbah/datasets/ukbb/ukbb_v3/var/chr1_22.afr_or_gbr.ldak.var.hapmap.cors.bim /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/prs/list_hg37_aou_top_mgbb.noukb.txt "ALT" "REF" /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/prs

## Step 1: estimate per-predictor heritabilities
${LDAK6} \
	--sum-hers ${per_snp_hers} \
	--summary ${gwas_summary_ldak_format} \
	--tagfile ${tagfile} \
	--matrix ${matrix} \
	--check-sums NO \
	--cutoff 0.01 \
	--max-threads ${maxThreads}

# Step 2:
# /medpop/esp2/mesbah/tools/LDAK5.2/ldak5.2.linux --mega-prs /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/bayesr_quickPRS.ldak.hg37.GWAMA.chr1_22.hasCH.MultiANC.AoU250k_TOPMed72k_MGBB53k.noukbb --summary /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/keep.ldak.hg37.GWAMA.chr1_22.hasCH.MultiANC.AoU250k_TOPMed72k_MGBB53k.noukbb.txt --ind-hers /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/quickPRS.ldak.hg37.GWAMA.chr1_22.hasCH.MultiANC.AoU250k_TOPMed72k_MGBB53k.noukbb.modified.ind.hers --cors /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/gbr.hapmap/gbr.hapmap --high-LD /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/gbr.hapmap/highld.snps --model bayesr --cv-proportion 0.1 --window-cm 1 --allow-ambiguous YES --extract /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/keep.ldak.hg37.GWAMA.chr1_22.hasCH.MultiANC.AoU250k_TOPMed72k_MGBB53k.noukbb.txt
##
${LDAK5} \
	--mega-prs ${mega_prs_weight} \
	--summary ${gwas_summary_ldak_format} \
	--ind-hers ${per_snp_hers}.ind.hers \
	--cors ${ref_cors} \
	--high-LD ${ref_highld} \
	--model bayesr \
	--cv-proportion 0.1 \
	--window-cm 1 \
	--extract ${gwas_summary_ldak_format} \
	--allow-ambiguous YES \
	--max-threads ${maxThreads}



