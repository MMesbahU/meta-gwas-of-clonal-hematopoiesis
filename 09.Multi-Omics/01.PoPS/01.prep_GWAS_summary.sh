#!/bin/bash

# for files in $(ls -l /medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/chr1_22.has*.21Aug2021_ukbEUR.tsv.gz | awk '{print $NF}'); do bash /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/PoPs/01.prep_GWAS_summary.sh ${files} /broad/hptmp/mesbah/ukb_chip/meta_gwas/pops/gwas/ukb/ukb_eur.magma_summary.$(basename ${files} ".gz"); done

input=${1}

output=${2}

echo -e "SNP\tP\tN" > ${output}
	# UKB
zcat ${input} | awk 'NR>1 && $13>=0.01 && $13<=0.99 && $24>=0.9{print $1"\t"$12"\t"($14+$18)}' >> ${output}

### from COJO input
## while read gene; do echo -e "SNP\tP\tN" > PoPs/PoPs_input/${gene}.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.tsv; done < <(echo -e "CHIP\nDNMT3A\nTET2")

## while read gene; do awk 'NR>1{print $1"\t"$7"\t"$8}' cojo/cojo_input/lifted_hg37.GWAMA.meta_${gene}.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.cojo.tsv >> PoPs/PoPs_input/${gene}.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.tsv; done < <(echo -e "CHIP\nDNMT3A\nTET2") &

## Meta-GWAS data
## while read files; do zcat ${files} | head -1 > /medpop/esp2/mesbah/projects/Meta_GWAS/sumstats/chr1_22.$(basename ${files} ".tsv.gz").tsv; for chr in {1..22}; do tabix ${files} ${chr} >> /medpop/esp2/mesbah/projects/Meta_GWAS/sumstats/chr1_22.$(basename ${files} ".tsv.gz").tsv; done; done < <(ls -l /medpop/esp2/mesbah/projects/Meta_GWAS/sumstats/lifted_hg37.GWAMA.meta_*.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.tsv.gz | awk '{print $NF}') &

## CHR	POS	REF	ALT	RSID	hg38MarkerID	EffectAllele	OtherAllele	EAF	BETA	SE	Pvalue	N	Direction	Het_P
while read files; do echo -e "SNP\tP\tN" > /broad/hptmp/mesbah/ukb_chip/PoPs/PoPs_input/chr1_22.$(basename ${files} ".tsv.gz").tsv; for chr in {1..22}; do tabix ${files} ${chr} | awk '{print $5"\t"$12"\t"$13}' >> /broad/hptmp/mesbah/ukb_chip/PoPs/PoPs_input/chr1_22.$(basename ${files} ".tsv.gz").tsv; done; done < <(ls -l /medpop/esp2/mesbah/projects/Meta_GWAS/sumstats/lifted_hg37.GWAMA.meta_*.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.tsv.gz | awk '{print $NF}') &



