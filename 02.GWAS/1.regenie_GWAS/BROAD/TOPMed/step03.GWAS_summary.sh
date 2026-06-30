#!/bin/bash

## while read ANC; do qsub -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/tmpdir -R y -l h_vmem=20G -l h_rt=20:00:00 -pe smp 1 -binding linear:1 -N step3.${ANC}.topmed74k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/BROAD/TOPMed/step03.GWAS_summary.sh ${ANC}; done < <(echo -e "MultiANC\nFemale\nMale\nAFR\nAMR\nEUR")


## combine chr1-22
# while read ANC; do while read traits; do echo -e "Name\tChr\tPos\tRef\tAlt\tTrait\tCohort\tModel\tEffect\tLCI_Effect\tUCI_Effect\tPval\tAAF\tNum_Cases\tCases_Ref\tCases_Het\tCases_Alt\tNum_Controls\tControls_Ref\tControls_Het\tControls_Alt\tInfo" | gzip -c > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/${ANC}/${ANC}.TOPMed74k.chr1_22.has${traits}.regenie.gz; while read lines; do zcat /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/${ANC}/step2/${ANC}.TOPMed74k.chr$(echo ${lines} | awk '{print $1}'):$(echo ${lines} | awk '{print $2}')_has${traits}.regenie.gz | awk 'NR>1'| gzip -c >> /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/${ANC}/${ANC}.TOPMed74k.chr1_22.has${traits}.regenie.gz; done </medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.nochr22_1_10163693.tsv; done < <(echo -e "CHvaf02\nCHvaf10\nDNMT3A\nTET2\nASXL1\nSF\nDDR"); done < <(echo -e "MultiANC\nFemale\nMale\nAFR\nAMR\nEUR") &
###
source /broad/software/scripts/useuse
##
use Tabix

use Bcftools

##
ANC=${1}
##
while read traits
do 
	echo -e "Name\tChr\tPos\tRef\tAlt\tTrait\tCohort\tModel\tEffect\tLCI_Effect\tUCI_Effect\tPval\tAAF\tNum_Cases\tCases_Ref\tCases_Het\tCases_Alt\tNum_Controls\tControls_Ref\tControls_Het\tControls_Alt\tInfo" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/${ANC}/${ANC}.TOPMed74k.chr1_22.has${traits}.regenie
	
	while read lines
	do 
		zcat /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/${ANC}/step2/${ANC}.TOPMed74k.chr$(echo ${lines} | awk '{print $1}'):$(echo ${lines} | awk '{print $2}')_has${traits}.regenie.gz | awk 'NR>1' >> /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/${ANC}/${ANC}.TOPMed74k.chr1_22.has${traits}.regenie 
	
	done </medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.nochr22_1_10163693.tsv

##  
	bgzip -f /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/${ANC}/${ANC}.TOPMed74k.chr1_22.has${traits}.regenie

## Index
	tabix -s 2 -b 3 -e 3 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/${ANC}/${ANC}.TOPMed74k.chr1_22.has${traits}.regenie.gz


done < <(echo -e "CHvaf02\nCHvaf10\nDNMT3A\nTET2\nASXL1\nSF\nDDR")

#####################################
##### Input for Meta-analyses  ######
#####################################

while read ANC; do while read phenos; do inDir="/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/${ANC}";outDir="/broad/hptmp/mesbah/dataset/ch_gwas/topmed/gwas/2024"; echo -e "SNPID\tRSID\tREF\tALT\tAAF\tBETA\tSE\tP\tN" | gzip -c > ${outDir}/GWAMA.${ANC}.TOPMed74k.chr1_22.has${phenos}.regenie.tsv.gz && zcat ${inDir}/${ANC}.TOPMed74k.chr1_22.has${phenos}.regenie.gz | sed '1d' | sed -e 's:;:\t:g' -e 's:REGENIE_BETA=::g' -e 's:REGENIE_SE=::g' | awk '{print "chr"$2":"$3":"$5":"$4"\t"$1"\t"$4"\t"$5"\t"$13"\t"$22"\t"$23"\t"$12"\t"($14+$18)}' | gzip -c >> ${outDir}/GWAMA.${ANC}.TOPMed74k.chr1_22.has${phenos}.regenie.tsv.gz; done< <(echo -e "CHvaf02\nCHvaf10\nDNMT3A\nTET2\nASXL1\nSF\nDDR"); done < <(echo -e "AFR\nAMR\nEUR\nFemale\nMale\nMultiANC") &








