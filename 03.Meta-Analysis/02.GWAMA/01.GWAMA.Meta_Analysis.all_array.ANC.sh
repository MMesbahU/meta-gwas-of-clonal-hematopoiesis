#!/bin/bash

source /broad/software/scripts/useuse

# use Tabix
################
## Make new dir
# mkdir -p /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/{Sex,anc}
#------
# prepare gwama input
# Female_Male
# while read Pheno; do echo -e "/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has${Pheno}.Female.ukbb450k_AoU250k_TOPMed72k_MGBB53k.out\n/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has${Pheno}.Male.ukbb450k_AoU250k_TOPMed72k_MGBB53k.out" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/ANC/${Pheno}.Female_Male.GWAMA.summary_list.txt;done < <(echo -e "CH\nCHvaf10\nDNMT3A\nTET2\nASXL1\nSF\nDDR")

# EUR_AFR
# while read Pheno; do echo -e "/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has${Pheno}.EUR.ukbb450k_AoU250k_TOPMed72k_MGBB53k.out\n/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has${Pheno}.AFR.ukbb450k_AoU250k_TOPMed72k_MGBB53k.out" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/ANC/${Pheno}.EUR_AFR.GWAMA.summary_list.txt;done < <(echo -e "CHvaf10\nASXL1\nSF\nDDR")
# while read Pheno; do echo -e "/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has${Pheno}.EUR.ukbb450k_AoU250k_TOPMed72k_MGBB53k_BioVU54k.out\n/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has${Pheno}.AFR.ukbb450k_AoU250k_TOPMed72k_MGBB53k.out" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/ANC/${Pheno}.EUR_AFR.GWAMA.summary_list.txt;done < <(echo -e "CH\nDNMT3A\nTET2")

# EUR + AFR + AMR
# while read Pheno; do echo -e "/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has${Pheno}.EUR.ukbb450k_AoU250k_TOPMed72k_MGBB53k.out\n/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has${Pheno}.AFR.ukbb450k_AoU250k_TOPMed72k_MGBB53k.out\n/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has${Pheno}.AMR.ukbb450k_AoU250k_TOPMed72k_MGBB53k.out" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/ANC/${Pheno}.EUR_AFR_AMR.GWAMA.summary_list.txt;done < <(echo -e "CHvaf10\nASXL1\nSF\nDDR")
# while read Pheno; do echo -e "/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has${Pheno}.EUR.ukbb450k_AoU250k_TOPMed72k_MGBB53k_BioVU54k.out\n/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has${Pheno}.AFR.ukbb450k_AoU250k_TOPMed72k_MGBB53k.out\n/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has${Pheno}.AMR.ukbb450k_AoU250k_TOPMed72k_MGBB53k.out" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/ANC/${Pheno}.EUR_AFR_AMR.GWAMA.summary_list.txt;done < <(echo -e "CH\nDNMT3A\nTET2")

# AFR + AMR
# while read Pheno; do echo -e "/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has${Pheno}.AFR.ukbb450k_AoU250k_TOPMed72k_MGBB53k.out\n/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has${Pheno}.AMR.ukbb450k_AoU250k_TOPMed72k_MGBB53k.out" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/ANC/${Pheno}.AFR_AMR.GWAMA.summary_list.txt;done < <(echo -e "CH\nDNMT3A\nTET2\nCHvaf10\nASXL1\nSF\nDDR")
#------------
# EUR + AMR
# while read Pheno; do echo -e "/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has${Pheno}.EUR.ukbb450k_AoU250k_TOPMed72k_MGBB53k.out\n/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has${Pheno}.AMR.ukbb450k_AoU250k_TOPMed72k_MGBB53k.out" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/ANC/${Pheno}.EUR_AMR.GWAMA.summary_list.txt;done < <(echo -e "CHvaf10\nASXL1\nSF\nDDR")
# while read Pheno; do echo -e "/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has${Pheno}.EUR.ukbb450k_AoU250k_TOPMed72k_MGBB53k_BioVU54k.out\n/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has${Pheno}.AMR.ukbb450k_AoU250k_TOPMed72k_MGBB53k.out" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/ANC/${Pheno}.EUR_AMR.GWAMA.summary_list.txt;done < <(echo -e "CH\nDNMT3A\nTET2")
#------

## get list of gwama summary
# ls /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/input/ANC/*summary_list.txt | awk '{print $NF}' > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/list_of_35_gwama_inputfiles.list

###### 
# mkdir -p /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/anc_stratified; 
# qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/tmpdir -t 1-$(wc -l /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/list_of_35_gwama_inputfiles.list |awk '{print $1}') -pe smp 4 -binding linear:4 -l h_vmem=30G -l h_rt=10:00:00 -N ANC_gwama_metaAnalyses /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/01.GWAMA.Meta_Analysis.all_array.ANC.sh rs_number n_samples reference_allele other_allele eaf beta se /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/list_of_35_gwama_inputfiles.list /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/anc_stratified /medpop/esp2/mesbah/tools/GWAMA_v2.2.2/GWAMA

###########################
# rs_number     reference_allele        other_allele    eaf     beta    se      beta_95L        beta_95U        p-value _-log10_p-value q_statistic     q_p-value       i2      n_studies       n_samples       effects

# chr1:10894:G:A        A       G       0.000722        0.470583        0.362322        -0.239569       1.180735        1.298797        0.193994        0.712212        2.795831        0.094510        0.642325        2       222525  +??+
# rs_number n_samples reference_allele other_allele eaf beta se
#---------------------------
name_marker=${1} # rs_number
name_n=${2} # n_samples
name_ea=${3} # reference_allele
name_nea=${4} # other_allele
name_eaf=${5} # eaf
name_beta=${6} # beta
name_se=${7} # se
#---------------------------
GWAMA_file_list=${8} # /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/list_of_96_gwama_inputfiles.list

outDir=${9}

## Run GWAMA
GWAMA=${10} # /medpop/esp2/mesbah/tools/GWAMA_v2.2.2/GWAMA
## input files
gwas_list=$(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${GWAMA_file_list} )

output_prefix=${outDir}/GWAMA.chr1_22.$(basename $(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${GWAMA_file_list} ) ".txt")
############################
# while read lines
# do 
#	gzip -d ${lines}.gz

# done <${gwas_list}
###################################################################

############### Task: Worked
## 1. Run Meta-Analysis
###########################

######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
##################################
# Ensure input file exists
if [[ ! -f "${gwas_list}" ]]; then
    echo "Error: Input file '${gwas_list}' not found."
    exit 1
fi

# Decompress all .gz files listed in ${gwas_list}
while read -r lines; do
    if ! gzip -d "${lines}.gz"; then
        echo "Error decompressing: ${lines}.gz"
        exit 1
    fi
    echo "Decompressed: ${lines}.gz"
done < "${gwas_list}"
#---------------------------------
# zcat /medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/chr1_22.has${traits}.21Aug2021_ukbEUR.tsv.gz | awk '(NR==1){print "SNPID\tN\tREF\tALT\tAAF\tBETA\tSE\tP"}(NR>1 && $2~/^[0-9]+/){print $1"\t"($14+$18)"\t"$4"\t"$5"\t"$13"\t"$22"\t"$23"\t"$12}' | awk '!seen[$1]++'| gzip -c > /medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/chr1_22.ukb200k.hg37_eur_${traits}.tsv.gz
##
## 

## Run GWAMA using GWAMA output
# GWAMA=/medpop/esp2/mesbah/tools/GWAMA_v2.2.2/GWAMA 
${GWAMA} \
	-i ${gwas_list} \
	-qt \
	--name_marker ${name_marker} \
	--name_n ${name_n} \
	--name_ea ${name_ea} \
	--name_nea ${name_nea} \
	--name_eaf ${name_eaf} \
	--name_beta ${name_beta} \
	--name_se ${name_se} \
	--indel_alleles \
	-o ${output_prefix}

## compress
#while read lines
#do
 #   gzip -f ${lines}

# done <${gwas_list}
# 
# gzip -f ${output_prefix}.out
# Recompress the files after processing
while read -r lines; do
	if ! gzip -f "${lines}"; then 
		echo "Error recompressing: ${lines}"
	        exit 1
	fi
	echo "Recompressed: ${lines}"
done < "${gwas_list}"

# Decompress the GWAMA output file
if ! gzip -f "${output_prefix}.out"; then
    echo "Error compressing: ${output_prefix}.out"
    exit 1
fi
echo "Compressed output file: ${output_prefix}.out.gz"
# 
###########


######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

