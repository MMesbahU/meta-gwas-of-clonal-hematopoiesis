#!/bin/bash

source /broad/software/scripts/useuse

### Get king software
# Get king v2.2.8  (released on May 10, 2022) : https://www.kingrelatedness.com/Download.shtml
# wget https://www.kingrelatedness.com/Linux-king.tar.gz -P /medpop/esp2/mesbah/tools
# tar -xzvf /medpop/esp2/mesbah/tools/Linux-king.tar.gz --directory /medpop/esp2/mesbah/tools/
##
#### king v2.2.8 needs following: 
# /lib64/libstdc++.so.6: version `CXXABI_1.3.8' 
# /lib64/libstdc++.so.6: version `CXXABI_1.3.9'

use GCC-5.2

############################# Clock Time: Start #####################################
echo -e "Job started at: $(date)"

Job_START=$(date +%s)
#######################################################################################


#######################  Run KING ########################################################
## Inputs
geno_bed_file=${1} # /broad/hptmp/mesbah/gwas/mgbb53k/qced.mgbb53k.autosomes.bed

out_prefix=${2} # /broad/hptmp/mesbah/gwas/mgbb53k/KIN3rd.qced.mgbb53k.autosomes

cpus=${3} # 8

##
## 3rd degree
## https://www.kingrelatedness.com/manual.shtml: "Please do not prune or filter any "good" SNPs that pass QC prior to any KING inference, unless the number of variants is too many to fit the computer memory, e.g., > 100,000,000 as in a WGS study, in which case rare variants can be filtered out. LD pruning is not recommended in KING."
##
/medpop/esp2/mesbah/tools/king \
	-b ${geno_bed_file} \
	--kinship \
	--degree 3 \
	--cpus ${cpus} \
	--prefix ${out_prefix}
#####################################################################################


############################### Clock Time: END ######################################
echo "Job ended at: $(date)" 

Job_END=$(date +%s)


echo $(( Job_END - Job_START)) | awk '{print "Total run time for KING: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'

#######################################################################################


