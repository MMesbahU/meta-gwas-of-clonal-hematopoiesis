#!/bin/bash

source /broad/software/scripts/useuse

use Bcftools

use Tabix

## qsub -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/gwas/aou/metadata -l h_rt=30:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 -R y -N combine_gwas_files /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/AoU/gwas_summary.SAS_EAS.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/gwas/aou
##

##
outDir=${1}

# Initialize the output file with the desired header
while read ANC; do
	
	while read pheno; do

		rm ${outDir}/${ANC}/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv{.gz,.gz.tbi}

		echo -e "SNPID\tCHROM\tPOS\tREF\tALT\tAAF\tBETA\tSE\tN\tP" > ${outDir}/${ANC}/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv

		# Loop through the 22 files corresponding to chr1 to chr22
		for file in ${outDir}/${ANC}/chr{1..22}.${ANC}.AoU_v71.has${pheno}.mac60.regenie.tsv.gz; do

		echo -e "${ANC}\t${pheno}: ${file}"
	        # Process each file to split SNPID and append to combined_output.tsv
        		zcat ${file} | awk 'NR > 1 { 
	        			   # Split SNPID into CHROM, POS, REF, ALT
		   	     	           split($1, snp, ":");
			        	   # Print SNPID, CHROM, POS, REF, ALT, followed by the remaining columns
				           print $1 "\t" snp[1] "\t" snp[2] "\t" snp[3] "\t" snp[4] "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6
					}' | sort -V -k1 >> ${outDir}/${ANC}/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv
   	       done

	       bgzip ${outDir}/${ANC}/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv

	       tabix -s 2 -b 3 -e 3 ${outDir}/${ANC}/chr1_22.${ANC}.AoU_v71.has${pheno}.regenie.tsv.gz

	done < <(echo -e "CH\nCHvaf10\nDNMT3A\nTET2\nASXL1\nSF\nDDR")

done < <(echo -e "EAS\nSAS")




