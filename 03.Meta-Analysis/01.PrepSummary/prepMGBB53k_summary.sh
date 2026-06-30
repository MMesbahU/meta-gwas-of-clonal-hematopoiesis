#!/bin/bash

###### 
# 53k
use Tabix

while read traits; do zcat /broad/hptmp/mesbah/gwas/mgbb53k/step2/MGBB53k.chr1\:1-24895642_hasCHIP.regenie.gz | head -1 > /medpop/esp/mesbah/GWAS_CHIP/MGBB53k/chr1_22.MGBB53k.${traits}.tsv; while read intervals; do zcat /broad/hptmp/mesbah/gwas/mgbb53k/step2/MGBB53k.chr$(echo ${intervals} | awk '{print $1}')\:$(echo ${intervals} | awk '{print $2}')_${traits}.regenie.gz | awk 'NR>1' >> /medpop/esp/mesbah/GWAS_CHIP/MGBB53k/chr1_22.MGBB53k.${traits}.tsv; done < <(cut -f1,2 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.tsv); done < <(awk '{print $1}' /broad/hptmp/mesbah/gwas/mgbb53k/NULL_MODEL_pred.list) &

## index 
for files in $(ls -lhv /medpop/esp/mesbah/GWAS_CHIP/MGBB53k/chr1_22.MGBB53k.*.tsv | awk '{print $NF}'); do bgzip -f ${files}; tabix -f -s 2 -b 3 -e 3 ${files}.gz; done &


## prepare for GWAMA: REF=Effect allele
# Name	Chr	Pos	Ref	Alt	Trait	Cohort	Model	Effect	LCI_Effect	UCI_Effect	Pval	AAF	Num_Cases	Cases_Ref	Cases_Het	Cases_Alt	Num_Controls	Controls_Ref	Controls_Het	Controls_Alt	Info
# chr1:10894:G:A	1	10894	A	G	hasCHIP	MGBB53k	ADD-WGR-FIRTH	0.593326	0.188408	1.86847	0.372449	0.999021	5550	0	8	5542	44755	0	41	44714	-0.522011	0.585283	INFO=0.345706	MAC=98.458824

while read traits; do zcat /medpop/esp/mesbah/GWAS_CHIP/MGBB53k/chr1_22.MGBB53k.has${traits}.tsv.gz | tr ';' '\t' | sed -e 's:REGENIE_BETA=::g' -e 's:REGENIE_SE=::g' | awk '(NR>1){SNPID=$1;REF=$4;ALT=$5;AAF=$13;BETA=$22;SE=$23;P=$12;N=($14+$18)}(NR==1){print "SNPID\tREF\tALT\tAAF\tBETA\tSE\tP\tN"}(NR>1){print SNPID"\t"REF"\t"ALT"\t"AAF"\t"BETA"\t"SE"\t"P"\t"N}' | awk '!seen[$1]++' > /broad/hptmp/mesbah/gwas/mgbb53k/mgbb53k_gwama.chr1_22.has${traits}.tsv; bgzip -f /broad/hptmp/mesbah/gwas/mgbb53k/mgbb53k_gwama.chr1_22.has${traits}.tsv; done < <(echo -e "CHIP\nDNMT3A\nTET2\nASXL1") &

## GWAMA input file list
while read traits; do echo -e "/medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/topmed2019.all_samples.has${traits}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/inputUKB/allUKB/lifted_hg38.chr1_22.has${traits}.11Aug2021_ukb200k_allsamples.INFO.tsv.gz\n/broad/hptmp/mesbah/gwas/mgbb53k/mgbb53k_gwama.chr1_22.has${traits}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_has${traits}_results_merged_subset.tsv.gz" > /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/01.PrepSummary/gwama_input.has${traits}.TOPMed2019UKB200kMGBB53kBioVU54k.txt; done < <(echo -e "CHIP\nDNMT3A\nTET2\nASXL1")



