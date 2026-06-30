#!/bin/bash


# ### !!! Note
#### REF allele is the effect allele in MGBB53k GWAS with --ref-first flag
#### To make ALT allele as EA, omit --ref-first flag
### !!!

######### Stratified GWAS
# zcat /medpop/esp/mesbah/GWAS_CHIP/MGBB53k/MultiAnc/step2/multi.MGBB53k.chr22\:10163694-20327387_hasCHvaf02.regenie.gz | head -2
## Name	Chr	Pos	Ref	Alt	Trait	Cohort	Model	Effect	LCI_Effect	UCI_Effect	Pval	AAF	Num_Cases	Cases_Ref	Cases_Het	Cases_Alt	Num_Controls	Controls_Ref	Controls_Het	Controls_Alt	Info
# chr22:10527916:T:G	22	10527916	T	G	hasCHvaf02	multi.MGBB53k	ADD-WGR-FIRTH	0.913158	0.675818	1.23385	0.554132	0.476895	5561	8	5553	0	44744	101	44643	0	REGENIE_BETA=-0.090847;REGENIE_SE=0.153567;INFO=0.871527;MAC=47980.439216;SCORE=-3.852264;SKATV=42.403929;LOG10P=0.256387

while read ANC; do while read phenos; do echo -e "SNPID\tREF\tALT\tAAF\tBETA\tSE\tP\tN\tINFO" | bgzip -c > /medpop/esp/mesbah/GWAS_CHIP/MGBB53k/has${phenos}.chr1_22.${ANC}.regenie.tsv.gz && for files in $(ls -lhv /medpop/esp/mesbah/GWAS_CHIP/MGBB53k/${ANC}/step2/${ANC}.MGBB53k.chr*:*_has${phenos}.regenie.gz | awk '{print $NF}'); do zcat ${files} | sed '1d' | sed -e 's:;:\t:g' -e 's:REGENIE_BETA=::g' -e 's:REGENIE_SE=::g' -e 's:INFO=::g' | awk '{print $1\t"$4"\t"$5"\t"$13"\t"$22"\t"$23"\t"$12"\t"($14+$18)"\t"$24}' | bgzip -c >> /medpop/esp/mesbah/GWAS_CHIP/MGBB53k/${ANC}/has${phenos}.chr1_22.${ANC}.regenie.tsv.gz; done; done < <(echo -e "CHvaf02\nCHvaf10\nDNMT3A\nTET2\nASXL1\nSF\nDDR");  done < <(echo -e "AMR\nEUR\nMale\nFemale") &

########

## EUR MGBB53k :
# while read trait; do qsub -wd /broad/hptmp/mesbah/gwas/mgbb53k/tmpdir -N mgbb53k_EUR.summary4GWAMA -l h_rt=10:00:00 -l h_vmem=10G -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/prepGWASsummary.MGB.sh /broad/hptmp/mesbah/gwas/mgbb53k/step2_EUR chr1_22.MGBB53k_EUR.has /broad/hptmp/mesbah/gwas/mgbb53k/step2_EUR/MGBB53k_EUR.chr ${trait}; done < <(echo -e "CHIP\nDNMT3A\nTET2")

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
# SNPID	REF	ALT	AAF	BETA	SE	P	N
############################################ MGBB
echo -e "SNPID\tChr\tPos_hg38\tREF\tALT\tTrait\tCohort\tModel\tEffect\tLCI_Effect\tUCI_Effect\tP\tAAF\tNum_Cases\tCases_Ref\tCases_Het\tCases_Alt\tNum_Controls\tControls_Ref\tControls_Het\tControls_Alt\tBETA\tSE\tINFO\tMAC\tN" | bgzip -c > ${outFile}

for files in $(ls -lhv ${inputPrefix}*_has${trait}.regenie.gz | awk '{print $NF}')

do 
	zcat ${files} | sed '1d' | sed -e 's:;:\t:g' -e 's:REGENIE_BETA=::g' -e 's:REGENIE_SE=::g' -e 's:INFO=::g' -e 's:MAC=::g' | awk '{print $0"\t"($14+$18)}' | bgzip -c >> ${outFile}

done
###########################################

