#!/bin/bash


## SNP effect: 
# /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.SNP_heritability/QuickPRS 1033 $ while read ANC; do ls /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/prs/mega/${ANC}/*.effects|awk '{print $NF}' | sort -V > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/prs/scores/list_megaprs_snp_effect.${ANC}.txt; done < <(echo -e "afr\ngbr")

# Run:
# while read ANC; do qsub -t 1-$(wc -l /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/prs/scores/list_megaprs_snp_effect.${ANC}.txt | awk '{print $1}') -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/tmpdir -N ${ANC}.ukb_scores -l h_rt=4:00:00 -l h_vmem=20G -pe smp 4 -binding linear:4 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.SNP_heritability/QuickPRS/get_scores.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/prs/scores/${ANC} /medpop/esp2/mesbah/datasets/ukbb/ukbb_v3/imp_v3/ukb_imp_chr1_22_v3 0 /medpop/esp2/mesbah/tools/LDAK6.0/ldak6.linux 4 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/sumHer/prs/scores/list_megaprs_snp_effect.${ANC}.txt; done < <(echo -e "afr\ngbr")

# /medpop/esp2/mesbah/tools/LDAK6.0/ldak6.linux --calc-scores /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/prs_score.hasCH.MultiANC.AoU250k_TOPMed72k_MGBB53k --scorefile /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/bayesr_quickPRS.ldak.hg37.GWAMA.chr1_22.hasCH.MultiANC.AoU250k_TOPMed72k_MGBB53k.noukbb.effects --bfile /medpop/esp2/mesbah/datasets/ukbb/ukbb_v3/ukb_imp_chr1_22_v3 --power 0 --max-threads 6

## https://dougspeed.com/profile-scores/
ldak_inpuDir_score=${1}

geno_file=${2} # /medpop/esp2/mesbah/datasets/ukbb/ukbb_v3/imp_v3/ukb_imp_chr1_22_v3

myPower=${3}  # 0

LDAK6=${4} # /medpop/esp2/mesbah/tools/LDAK6.0/ldak6.linux | /medpop/esp2/mesbah/tools/LDAK5.2/ldak5.2.linux

maxThreads=${5}

list_mega_prs_weight=${6}

####
### get GWAS file from t array file list
mega_prs_weight=$(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${list_mega_prs_weight} )

score_outfile=${ldak_inpuDir_score}/prs_score.$(basename $(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${list_mega_prs_weight} ) ".effects")
####


##
${LDAK6} \
	--calc-scores ${score_outfile} \
	--scorefile ${mega_prs_weight} \
	--bfile ${geno_file} \
	--power ${myPower} \
	--max-threads ${maxThreads}

