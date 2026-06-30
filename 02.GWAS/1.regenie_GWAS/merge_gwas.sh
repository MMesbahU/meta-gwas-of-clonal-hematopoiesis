#!/bin/bash

## qsub -wd /medpop/esp/mesbah/GWAS_CHIP/MGBB53k/tmpdir -R y -l h_rt=20:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 -N combine_summary /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/merge_gwas.sh



## rename files: for files in $(ls multi.MGBB53k.chr* | awk '{print $NF}'); do mv ${files} $(echo $files | sed 's:multi:MultiAnc:g'); done
## for files in $(ls Female_XX.MGBB53k.chr* | awk '{print $NF}'); do mv ${files} $(echo $files | sed 's:Female_XX:Female:g'); done


source /broad/software/scripts/useuse

use Tabix


##
while read ANC; do 
	while read phenos; do
	echo -e "SNPID\tREF\tALT\tAAF\tBETA\tSE\tP\tN\tINFO" | bgzip -c > /medpop/esp/mesbah/GWAS_CHIP/MGBB53k/${ANC}/has${phenos}.chr1_22.${ANC}.regenie.tsv.gz
	for files in $(ls -lhv /medpop/esp/mesbah/GWAS_CHIP/MGBB53k/${ANC}/step2/${ANC}.MGBB53k.chr*:*_has${phenos}.regenie.gz | awk '{print $NF}'); do
	zcat ${files} | sed '1d' | sed -e 's:;:\t:g' -e 's:REGENIE_BETA=::g' -e 's:REGENIE_SE=::g' -e 's:INFO=::g' | awk '{print $1"\t"$4"\t"$5"\t"$13"\t"$22"\t"$23"\t"$12"\t"($14+$18)"\t"$24}' | bgzip -c >> /medpop/esp/mesbah/GWAS_CHIP/MGBB53k/${ANC}/has${phenos}.chr1_22.${ANC}.regenie.tsv.gz
	done; done < <(echo -e "CHvaf02\nCHvaf10\nDNMT3A\nTET2\nASXL1\nSF\nDDR"); done < <(echo -e "AFR\nSAS\nEAS\nMale\nFemale\nMultiAnc")


#### echo -e "AMR\nEUR"

