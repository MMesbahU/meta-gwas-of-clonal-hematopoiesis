#!/bin/bash


### Oct 2024

##


## WB: 
## while read trait; do qsub -wd /medpop/esp/mesbah/GWAS_CHIP/inputUKB/ukb450k/250k/tmpdir -N ukb250k_WB.sum4GWAMA -l h_rt=10:00:00 -l h_vmem=10G -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/prepGWASsummary.UKB.sh /medpop/esp/mesbah/GWAS_CHIP/inputUKB/ukb450k/250k/hg37 chr1_22.ukb250k_WB.has /broad/hptmp/mesbah/gwas/ukb450k/250k.step2/ukb250k_WB/UKB250k_WB.chr ${trait}; done < <(echo -e "CHIP\nDNMT3A\nTET2")

## multi-ancestry:
## while read trait; do qsub -wd /medpop/esp/mesbah/GWAS_CHIP/inputUKB/ukb450k/250k/tmpdir -N ukb250k_multi_ancestry.sum4GWAMA -l h_rt=10:00:00 -l h_vmem=10G -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/prepGWASsummary.UKB.sh /medpop/esp/mesbah/GWAS_CHIP/inputUKB/ukb450k/250k/hg37 chr1_22.ukb250k_multiAncestry.has /broad/hptmp/mesbah/gwas/ukb450k/250k.step2/additive/UKB250k.chr ${trait}; done < <(echo -e "CHIP\nDNMT3A\nTET2")
##
source /broad/software/scripts/useuse
use Tabix
##
outDir=${1}
out_prefix=${2}
inputPrefix=${3}
trait=${4}
outFile=${outDir}/${out_prefix}${trait}.tsv.gz # chr1_22.ukb450.has
# echo -e "Name\tChr\tPos_hg19\tRef\tAlt\tTrait\tCohort\tModel\tEffect\tLCI_Effect\tUCI_Effect\tPval\tAAF\tNum_Cases\tCases_Ref\tCases_Het\tCases_Alt\tNum_Controls\tControls_Ref\tControls_Het\tControls_Alt\tREGENIE_BETA\tREGENIE_SE\tINFO\tMAC" | gzip -c > chr1_22.additive_hasCHIP.regenie.gz
# gsutil cat gs://ukbb_v2/projects/muddin/chip_gwas/UKB200k_11Aug2021/additive/Chr1_additive_hasCHIP.regenie.gz| zcat | sed '1d' | sed -e 's:;:\t:g' -e 's:REGENIE_BETA=::g' -e 's:REGENIE_SE=::g' -e 's:INFO=::g' -e 's:MAC=::g' | gzip -c > chr1_22.additive_hasCHIP.regenie.gz

###################################################### 
## 450k
# while read gene_name 
# do
#	echo -e "Name\tChr\tPos_hg19\tRef\tAlt\tTrait\tCohort\tModel\tEffect\tLCI_Effect\tUCI_Effect\tPval\tAAF\tNum_Cases\tCases_Ref\tCases_Het\tCases_Alt\tNum_Controls\tControls_Ref\tControls_Het\tControls_Alt\tREGENIE_BETA\tREGENIE_SE\tINFO\tMAC" | bgzip -c > /broad/hptmp/mesbah/gwas/ukb450k/450k.step2/chr1_22.has${gene_name}.regenie.gz
#
#	for files in $(ls -lhv /broad/hptmp/mesbah/gwas/ukb450k/450k.step2/additive/UKB450k.chr*_has${gene_name}.regenie.gz | awk '{print $NF}')
#	do
#		zcat ${files} | sed '1d' | sed -e 's:;:\t:g' -e 's:REGENIE_BETA=::g' -e 's:REGENIE_SE=::g' -e 's:INFO=::g' -e 's:MAC=::g' | bgzip -c >> /broad/hptmp/mesbah/gwas/ukb450k/450k.step2/chr1_22.has${gene_name}.regenie.gz
#	done
#
# done < <(echo -e "CHIP\nDNMT3A\nTET2\nASXL1")
##
###################################################### 

############################################ 250k
echo -e "Name\tChr\tPos_hg19\tREF\tALT\tTrait\tCohort\tModel\tEffect\tLCI_Effect\tUCI_Effect\tP\tAAF\tNum_Cases\tCases_Ref\tCases_Het\tCases_Alt\tNum_Controls\tControls_Ref\tControls_Het\tControls_Alt\tBETA\tSE\tINFO\tMAC\tN" | bgzip -c > ${outFile}

for files in $(ls -lhv ${inputPrefix}*_has${trait}.regenie.gz | awk '{print $NF}')

do 
	zcat ${files} | sed '1d' | sed -e 's:;:\t:g' -e 's:REGENIE_BETA=::g' -e 's:REGENIE_SE=::g' -e 's:INFO=::g' -e 's:MAC=::g' | awk '{print $0"\t"($14+$18)}' | bgzip -c >> ${outFile}

done
###########################################

