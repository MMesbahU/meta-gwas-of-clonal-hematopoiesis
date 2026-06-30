#!/bin/bash

## qsub -wd /broad/hptmp/mesbah/dataset/ch_gwas/ukbb/MultiANC -b y -N cp_gcp -l h_rt=30:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 -S /bin/bash 'source /broad/software/scripts/useuse && use Google-Cloud-SDK && cd /broad/hptmp/mesbah/dataset/ch_gwas/ukbb/MultiANC && gsutil -m cp -r UKB450k/MultiANC gs://muddin/chip_gwas/'

### Oct 2024
# NOTE: have to add "chr"
# mkdir -p /broad/hptmp/mesbah/dataset/ch_gwas/ukbb/MultiANC/tmpdir; while read lines; do qsub -R y -t 1-8 -wd /broad/hptmp/mesbah/dataset/ch_gwas/ukbb/MultiANC/tmpdir -N prep_sum.${lines} -l h_rt=30:00:00 -l h_vmem=10G -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/RAP/ukb450k.gwas_summary.prep.sh /broad/hptmp/mesbah/dataset/ch_gwas/ukbb/MultiANC/combined 'chr1_22.MultiANC' /broad/hptmp/mesbah/dataset/ch_gwas/ukbb/MultiANC ${lines} 'MultiAnc'; done < <(echo -e "ukb200k_N193342\nukb250k_N243350")
##
source /broad/software/scripts/useuse

##
outDir=${1}
out_prefix=${2} # 'chr1_22.MultiANC
inputDir=${3} # /broad/hptmp/mesbah/dataset/ch_gwas/ukbb/MultiANC/ 
# ukb250k_N243350/chr22_MultiAnc.ukb250k_N243350.22_20327388-30491080_hasCHvaf10.regenie.gz
cohorts=${4} # ukb250k_N243350
ANC=${5} # MultiAnc
####
trait_num=${SGE_TASK_ID} # 1-8

trait=$(echo -e 'SF\nDDR\nASXL1\nTET2\nDNMT3A\nCHvaf10\nCHvaf05\nCH' | awk -v mynum=${trait_num} 'NR==mynum{print $1}')
# 
outFile=${outDir}/${out_prefix}.${cohorts}.has${trait}.regenie.tsv

###################################################### 
############################################ 

# File Header
echo -e "SNPID\tREF\tALT\tAAF\tBETA\tSE\tP\tN\tINFO" > ${outFile}
##

echo -e "${cohorts}\t${trait}\n"

echo -e "Processing file:\n"
#
for files in $(ls -lh ${inputDir}/${cohorts}/chr*_${ANC}.${cohorts}.*_has${trait}.regenie.gz | awk '{print $NF}' | sort -V)

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

