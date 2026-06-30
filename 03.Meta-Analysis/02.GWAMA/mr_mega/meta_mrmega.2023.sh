#!/bin/bash

######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
#########

## Worked: /home/jupyter/workspaces/detectionofclonalhematopoiesisofindeterminatepotentialchip/MR_MEGA_v0.2/MR-MEGA -i mgb_aou_ch.txt -o mrmega.meta.mgb_aou_ch --qt --pc 1 --no_std_names --name_n N --name_se SE --name_beta Beta --name_eaf EAF --name_ea EA --name_nea NEA --name_marker SNPID --name_pos POS --name_chr CHR

MRMEGA=${1}
summary_list=${2}
nPC=${3}
outPrefix=${4}

###
${MRMEGA} \
	-i ${summary_list} \
	-o ${outPrefix} \
	--qt \
	--pc ${nPC} \
	--no_std_names \
	--name_n N \
	--name_se SE \
	--name_beta Beta \
	--name_eaf EAF \
	--name_ea EA \
	--name_nea NEA \
	--name_marker SNPID \
	--name_pos POS \
	--name_chr CHR


## 
gzip ${outPrefix}.result

######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########


