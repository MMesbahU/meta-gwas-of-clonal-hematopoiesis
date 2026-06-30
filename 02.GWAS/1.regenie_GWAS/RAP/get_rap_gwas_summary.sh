#!/bin/bash

## qsub -wd /broad/hptmp/mesbah/dataset/ch_gwas/ukbb/MultiANC -N get_summary -R y -l h_vmem=10G -l h_rt=20:00:00 -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/RAP/get_rap_gwas_summary.sh

### To load dx; we need python3
source /broad/software/scripts/useuse

use Anaconda3

## /home/unix/muddin/.local/bin/dx login 
# provide user name
# provide pass
# /home/unix/muddin/.local/bin/dx ls ukbb_CH/2024/gwas/step2/MultiAnc/{ukb200k_N193342,ukb250k_N243350}


### September 2024
while read cohorts
do 
	/home/unix/muddin/.local/bin/dx download -a -f --recursive ukbb_CH/2024/gwas/step2/MultiAnc/${cohorts} -o /broad/hptmp/mesbah/dataset/ch_gwas/ukbb/MultiANC/ 
done < <(echo -e "ukb200k_N193342\nukb250k_N243350")


##

########################### March 2024
######## Stratified GWAS Summary statis in UKB
## while read trait; do /home/unix/muddin/.local/bin/dx download -a -f ukbb_CH/2024/gwas/step2/${trait}/*.gz -o /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/step2/${trait}/; done < <(echo -e "AFR\nAMR\nSAS\nEAS") 1>>/medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/step2/rap.03_07_2024.download.log 2>>/medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/step2/rap.03_07_2024.download.err &
# 
# while read trait; do /home/unix/muddin/.local/bin/dx download --recursive ukbb_CH/2024/gwas/step2/${trait}; done < <(echo -e "Male\nFemale\nEUR\nAFR\nAMR\nSAS\nEAS") 1>>/medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/step2/rap.03_07_2024.download.log 2>>/medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/step2/rap.03_07_2024.download.err &

### Summary Stats: 
## note: need to add "chr" in ukbb GWAS hg38 mapping
## # while read ANC; do while read phenos; do echo -e "SNPID\tREF\tALT\tAAF\tBETA\tSE\tP\tN\tINFO" | bgzip -c > /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/step2/summary/has${phenos}.chr1_22.${ANC}.regenie.tsv.gz && for files in $(ls -lhv /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/step2/${ANC}/chr*_${ANC}.*_has${phenos}.regenie.gz | awk '{print $NF}'); do zcat ${files} | sed '1d' | sed -e 's:;:\t:g' -e 's:REGENIE_BETA=::g' -e 's:REGENIE_SE=::g' -e 's:INFO=::g' -e 's:MAC=::g' | awk '{print "chr"$2":"$3":"$4":"$5"\t"$4"\t"$5"\t"$13"\t"$22"\t"$23"\t"$12"\t"($14+$18)"\t"$24}' | bgzip -c >> /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/step2/summary/has${phenos}.chr1_22.${ANC}.regenie.tsv.gz; done; done < <(echo -e "CH\nCHvaf10\nDNMT3A\nTET2\nASXL1\nSF\nDDR");  done < <(echo -e "AFR\nAMR\nSAS\nEAS\nEUR\nMale\nFemale")

############################

## ukb200k
# for chr in {1..22}; do /home/unix/muddin/.local/bin/dx download ukbb_CH/ch_gwas_2023/step2/chr${chr}_ukb200k.*.regenie.gz; done


## ukb250k
# for chr in {1..22}; do /home/unix/muddin/.local/bin/dx download ukbb_CH/ch_gwas_2023/step2/chr${chr}_ukb250k.*.regenie.gz; done


##
# zcat ukb200k/chr22_ukb200k.22_10163694-20327387_hasCH.regenie.gz  | sed '1d' | sed -e 's:;:\t:g' -e 's:REGENIE_BETA=::g' -e 's:REGENIE_SE=::g' -e 's:INFO=::g' -e 's:MAC=::g' | awk '{print $2":"$3":"$4":"$5"\t"$4"\t"$5"\t"$13"\t"$22"\t"$23"\t"($14+$18)"\t"$12"\t"$24}'

# use Tabix
	## UKB200k
# while read phenos; do echo -e "SNPID\tREF\tALT\tAAF\tBETA\tSE\tP\tN\tINFO" | bgzip -c > /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2023_ukbb/rap_gwas/has${phenos}.chr1_22.ukb200k.regenie.tsv.gz && for files in $(ls -lhv /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2023_ukbb/rap_gwas/ukb200k/chr*_ukb200k.*_has${phenos}.regenie.gz | awk '{print $NF}'); do zcat ${files} | sed '1d' | sed -e 's:;:\t:g' -e 's:REGENIE_BETA=::g' -e 's:REGENIE_SE=::g' -e 's:INFO=::g' -e 's:MAC=::g' | awk '{print $2":"$3":"$4":"$5"\t"$4"\t"$5"\t"$13"\t"$22"\t"$23"\t"$12"\t"($14+$18)"\t"$24}' | bgzip -c >> /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2023_ukbb/rap_gwas/has${phenos}.chr1_22.ukb200k.regenie.tsv.gz; done; done< <(echo -e "CH\nDNMT3A\nTET2\nDTA") &


	## UKB250k
# while read phenos; do echo -e "SNPID\tREF\tALT\tAAF\tBETA\tSE\tP\tN\tINFO" | bgzip -c > /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2023_ukbb/rap_gwas/has${phenos}.chr1_22.ukb250k.regenie.tsv.gz && for files in $(ls -lhv /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2023_ukbb/rap_gwas/ukb250k/chr*_ukb250k.*_has${phenos}.regenie.gz | awk '{print $NF}'); do zcat ${files} | sed '1d' | sed -e 's:;:\t:g' -e 's:REGENIE_BETA=::g' -e 's:REGENIE_SE=::g' -e 's:INFO=::g' -e 's:MAC=::g' | awk '{print $2":"$3":"$4":"$5"\t"$4"\t"$5"\t"$13"\t"$22"\t"$23"\t"$12"\t"($14+$18)"\t"$24}' | bgzip -c >> /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2023_ukbb/rap_gwas/has${phenos}.chr1_22.ukb250k.regenie.tsv.gz; done; done< <(echo -e "CH\nDNMT3A\nTET2\nDTA") &



## add "chr"
# while read phenos; do echo -e "SNPID\tREF\tALT\tAAF\tBETA\tSE\tP\tN\tINFO" > /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2023_ukbb/rap_gwas/has${phenos}.chr1_22.ukb200k.regenie.tsv && zcat /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2023_ukbb/rap_gwas/has${phenos}.chr1_22.ukb200k.regenie.tsv.gz | awk 'NR>1{print "chr"$0}' >> /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2023_ukbb/rap_gwas/has${phenos}.chr1_22.ukb200k.regenie.tsv; done< <(echo -e "CH\nDNMT3A\nTET2\nDTA") &

# while read pheno; do echo -e "SNPID\tREF\tALT\tAAF\tBETA\tSE\tP\tN\tINFO" > /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2023_ukbb/rap_gwas/has${pheno}.chr1_22.ukb250k.regenie.tsv && zcat /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2023_ukbb/rap_gwas/has${pheno}.chr1_22.ukb250k.regenie.tsv.gz | awk 'NR>1{print "chr"$0}' >> /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2023_ukbb/rap_gwas/has${pheno}.chr1_22.ukb250k.regenie.tsv; done< <(echo -e "CH\nDNMT3A\nTET2\nDTA") &

# gzip -f /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2023_ukbb/rap_gwas/has${phenos}.chr1_22.ukb200k.regenie.tsv

