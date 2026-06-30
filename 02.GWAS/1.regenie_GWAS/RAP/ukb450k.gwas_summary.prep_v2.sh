#!/bin/bash


### Oct 2024
# NOTE: have to add "chr"
# qsub -R y -wd /broad/hptmp/mesbah/dataset/ch_gwas/tmpdir -N prep_sum.Female_DDR -l h_rt=10:00:00 -l h_vmem=10G -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/RAP/ukb450k.gwas_summary.prep_v2.sh /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb 'UKB450k.chr1_22' /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/step2 DDR Female
##
source /broad/software/scripts/useuse

##
outDir=${1} # /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/input/ukbb

out_prefix=${2} # UKB450k.chr1_22

inputDir=${3} # /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/step2

trait=${4} # DDR

ANC=${5} # Fenale
####
# trait_num=${SGE_TASK_ID} # 1-8

# trait=$(echo -e 'SF\nDDR\nASXL1\nTET2\nDNMT3A\nCHvaf10\nCHvaf05\nCH' | awk -v mynum=${trait_num} 'NR==mynum{print $1}')
# 
outFile=${outDir}/has${trait}.${out_prefix}.${ANC}.regenie.tsv

# /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/step2/Female/chr22_Female.22_40654775-50818468_hasDDR.regenie.gz head hasDDR.UKB450k.chr1_22.Female.regenie.tsv

###################################################### 
############################################ 

# File Header
echo -e "SNPID\tREF\tALT\tAAF\tBETA\tSE\tP\tN\tINFO" > ${outFile}
##

echo -e "${ANC}\t${trait}\n"

echo -e "Processing file:\n"
#
for files in $(ls -lh ${inputDir}/${ANC}/chr*_${ANC}.*_has${trait}.regenie.gz | awk '{print $NF}' | sort -V)

do 
	echo "${files}"

	zcat ${files} | sed '1d' | cut -f2-5,12,13,14,18,22 | \
	awk -F'\t' '{
	    # Create SNPID by merging columns 1, 2, 3, and 4 with ":"
	    # add "chr"
	    snpid = "chr"$1":"$2":"$3":"$4;

	    # Sum columns 7 and 8 to create N
	    N = $7 + $8;

	    # Extract beta, se, and info from column 9
	    split($9, fields, ";");

	    # Initialize variables for REGENIE_BETA, REGENIE_SE, and INFO
	    beta=""; se=""; info="";

	    # Loop through each field and extract the values for REGENIE_BETA, REGENIE_SE, and INFO
	    for(i=1;i<=length(fields);i++) {
	        if (fields[i] ~ /REGENIE_BETA=/) {
	            split(fields[i], a, "="); beta=a[2]
	        }
	        if (fields[i] ~ /REGENIE_SE=/) {
	            split(fields[i], a, "="); se=a[2]
	        }
	        if (fields[i] ~ /INFO=/) {
	            split(fields[i], a, "="); info=a[2]
	        }
																     }

																     # Print the final output in the desired order:
																     # SNPID, REF, ALT, AAF, beta, se, P, N, info
																     printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", snpid, $3, $4, $6, beta, se, $5, N, info
								
      }' >> ${outFile}

done

###########################################

