#!/bin/bash

# qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas -N gcp_TOPgwas -l h_rt=2:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/BROAD/TOPMed/step05.cp_topmed_gwas_summary.sh "/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas" "gs://muddin/chip_gwas"
#
source /broad/software/scripts/useuse

use Google-Cloud-SDK

inDir=${1} # /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/
outDir=${2} # gs://muddin/chip_gwas
## 
while read ANC
do 
	while read pheno 
	
	do 
		gsutil -m cp ${inDir}/${ANC}/${ANC}.TOPMed74k.chr1_22.has${pheno}.regenie.gz ${outDir}/
	
	done < <(echo -e "ASXL1\nCHvaf10\nSF\nDDR")

done < <(echo -e "MultiANC\nFemale\nMale\nEUR\nAFR\nAMR" )


