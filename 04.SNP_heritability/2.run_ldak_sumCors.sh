#!/bin/bash


# /medpop/esp2/mesbah/tools/LDAK5.2/ldak5.2.linux --summary /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/in_ldak_eaf001_min2Studies.8March2022.lifted_hg37.GWAMA.meta_CHIP.allsamples_topmed_ukbb_mgbb_BioVU.txt --summary2 /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/input_MPN.GWAS.txt --tagfile /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/ldak.thin.genotyped.gbr.tagging --sum-cors /broad/hptmp/mesbah/ukb_chip/h2/geneticCor/genCor.CHIP_MPN --check-sums NO

# MPN downloded from GWAScatalog: /medpop/esp2/mesbah/tools/LDAK5.2/ldak5.2.linux --summary /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/in_ldak_eaf001_min2Studies.8March2022.lifted_hg37.GWAMA.meta_CHIP.allsamples_topmed_ukbb_mgbb_BioVU.txt --summary2 /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/input_MPN_GWAS.GCST90000032_GRCh37.txt --tagfile /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/ldak.thin.genotyped.gbr.tagging --sum-cors /broad/hptmp/mesbah/ukb_chip/h2/geneticCor/genCor.CHIP_MPN.GCST90000032_GRCh37 --check-sums NO 1>>/broad/hptmp/mesbah/ukb_chip/h2/geneticCor/chip_mpnGCST90000032_GRCh37.log 2>> /broad/hptmp/mesbah/ukb_chip/h2/geneticCor/chip_mpnGCST90000032_GRCh37.err &

## LTL GWAS: /medpop/esp2/mesbah/tools/LDAK5.2/ldak5.2.linux --summary /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/in_ldak_eaf001_min2Studies.8March2022.lifted_hg37.GWAMA.meta_CHIP.allsamples_topmed_ukbb_mgbb_BioVU.txt --summary2 /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/input_UKB_telomere_gwas_summarystats.N78592_GRCh37.txt --tagfile /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/ldak.thin.genotyped.gbr.tagging --sum-cors /broad/hptmp/mesbah/ukb_chip/h2/geneticCor/genCor.CHIP_UKB_telomereN78592 --check-sums NO 1>>/broad/hptmp/mesbah/ukb_chip/h2/geneticCor/chip_UKB_telomereN78592.log 2>>/broad/hptmp/mesbah/ukb_chip/h2/geneticCor/chip_UKB_telomereN78592.errr

#### http://dougspeed.com/pre-computed-tagging-files/
# Get precomputed Tagging SNPs
# LDAK-Thin Model: wget -P /broad/hptmp/mesbah/ukb_chip/h2/ https://genetics.ghpc.au.dk/doug/ldak.thin.genotyped.gbr.tagging.gz
# BLD-LDAK-Lite-Alpha Model: wget https://genetics.ghpc.au.dk/doug/bld.ldak.lite.alpha.hapmap.gbr.tagging.gz
# BLD-LDAK Model: wget -P /broad/hptmp/mesbah/ukb_chip/ https://genetics.ghpc.au.dk/doug/bld.ldak.genotyped.gbr.tagging.gz
# 

## Sumher: http://dougspeed.com/snp-heritability/
sum_trait1=${1} # to specify the file containing the summary statistics.
sum_trait2=${2}
taggFile=${3}  # to specify the tagging file
out_file_prefix=${4}
### If the summary statistics come from analysing a binary phenotype, then you can use --prevalence <float> and --ascertainment <float> to specify the proportion of cases in the population and in the GWAS; LDAK will then also report estimates of variance explained on the liability scale.

## run ldak
/medpop/esp2/mesbah/tools/LDAK5.2/ldak5.2.linux \
	--summary ${sum_trait1} \
	--summary2 ${sum_trait2}
	--tagfile ${taggFile} \
	--sum-cors ${out_file_prefix} \
	--check-sums NO


### Correlation
## while read refs; do /medpop/esp2/mesbah/tools/LDAK5.2/ldak5.2.linux --summary /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/in_ldak_eaf001_min2Studies.8March2022.lifted_hg37.GWAMA.meta_CHIP.allsamples_topmed_ukbb_mgbb_BioVU.txt --summary2 /broad/hptmp/mesbah/ukb_chip/ltl_gwas/in_ldak.UKB_telomere_gwas_summarystats.NatGen_N472174.GRCh37.txt --tagfile /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/ldak-thin/ldak.thin.genotyped.${refs}.tagging --sum-cors /broad/hptmp/mesbah/ukb_chip/h2/genCor/m4.CHIP.m4All.LTL_ukbN472174.ldak_thin.geno_${refs} --check-sums NO 1>>/broad/hptmp/mesbah/ukb_chip/h2/genCor/cor.chip_ltl.log 2>>/broad/hptmp/mesbah/ukb_chip/h2/genCor/cor.chip_ltl.err; done < <(echo -e "gbr\nsas\neas\nafr") &
## Eur only: while read refs; do /medpop/esp2/mesbah/tools/LDAK5.2/ldak5.2.linux --summary /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/in_ldak_eaf001_min2Studies.8March2022.lifted_hg37.GWAMA.EUR_only_CHIP_GWAS.TopMedUKbbMGbbBioVU.txt --summary2 /broad/hptmp/mesbah/ukb_chip/ltl_gwas/in_ldak.UKB_telomere_gwas_summarystats.NatGen_N472174.GRCh37.txt --tagfile /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/ldak-thin/ldak.thin.genotyped.${refs}.tagging --sum-cors /broad/hptmp/mesbah/ukb_chip/h2/genCor/CHIP.m4_eurOnly.LTL_ukbN472174.ldak_thin.geno_${refs} --check-sums NO 1>>/broad/hptmp/mesbah/ukb_chip/h2/genCor/cor.chip_ltl.log 2>>/broad/hptmp/mesbah/ukb_chip/h2/genCor/cor.chip_ltl.err; done < <(echo -e "gbr\nsas\neas\nafr") &

