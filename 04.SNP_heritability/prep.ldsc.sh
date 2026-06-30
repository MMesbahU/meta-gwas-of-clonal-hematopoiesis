#!/bin/bash


## qsub -wd /broad/hptmp/mesbah/ukb_chip/tmpdir -R y -l h_rt=20:00:00 -l h_vmem=20G -N prep_ldsc /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/09.Heritability/prep.ldsc.sh

## /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/09.Heritability/run.ldsc.sh /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/summary_for_ldsc.topmed.hg37.CHIP.EUR.results.tsv /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/eur_w_ld_chr /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/eur_w_ld_chr /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/Her_ldsc.topmed.hg37.CHIP.EUR.results

source /broad/software/scripts/useuse

use Anaconda

# cd /medpop/esp2/mesbah/tools/ldsc

# conda env create --file environment.yml
source activate ldsc
# source deactivate

# https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation
# gwas_summary=${1} # /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/summary_for_ldsc.topmed.hg37.CHIP.EUR.results.tsv

# ref_ld_dir={2} # /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/eur_w_ld_chr

# w_ld_dir=${3} # /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/eur_w_ld_chr

# outHer=${4} # 
	## get LD for EUR 
	# wget -P /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/ https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
	# tar -jxvf /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/eur_w_ld_chr.tar.bz2

## data formating
# /medpop/esp2/mesbah/tools/ldsc/munge_sumstats.py --out ldsc.topmed.hg37.CHIP.EUR.results --sumstats /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/summary_for_ldsc.topmed.hg37.CHIP.EUR.results.tsv 

	## Run h2 in EUR
#/medpop/esp2/mesbah/tools/ldsc/ldsc.py \
#	--h2 ${gwas_summary} \
#	--ref-ld-chr ${ref_ld_dir}/ \
#	--w-ld-chr ${w_ld_dir}/ \
#	--out ${outHer}

## /medpop/esp2/mesbah/tools/ldsc/ldsc.py --h2 /broad/hptmp/mesbah/ukb_chip/h2/ldsc.topmed.hg37.CHIP.EUR.results.sumstats.gz --ref-ld-chr /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/eur_w_ld_chr/ --w-ld-chr /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/eur_w_ld_chr/ --out /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/Her_ldsc.topmed.hg37.CHIP.EUR.results
## liability scale
## /medpop/esp2/mesbah/tools/ldsc/ldsc.py --h2 /broad/hptmp/mesbah/ukb_chip/h2/ldsc.topmed.hg37.CHIP.EUR.results.sumstats.gz --ref-ld-chr /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/eur_w_ld_chr/ --w-ld-chr /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/eur_w_ld_chr/ --out /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/Her_ldsc.topmed.hg37.CHIP.EUR.results --samp-prev 0.05 --pop-prev 0.05


### Prepare UKB CHIP GWAS
# /medpop/esp2/mesbah/tools/ldsc/munge_sumstats.py --sumstats /medpop/esp2/mesbah/projects/Meta_GWAS/Sep2021/Zhi_allUKB/chr1_22.hasCHIP.11Aug2021_ukb200k.allUKB_MAF001_INFO_03.tsv.gz --out /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/ldsc.chr1_22.hasCHIP.11Aug2021_ukb200k.allUKB_MAF001_INFO_03 --a1 EffectAllele --a2 OtherAllele --p Pvalue --N-col N --snp RSID --frq EAF

## Prep LTL GWAS: 
/medpop/esp2/mesbah/tools/ldsc/munge_sumstats.py \
	--sumstats /medpop/esp/mesbah/GWAS_CHIP/ltl_gwas/UKB_telomere_gwas_summarystats.tsv.gz \
	--out /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/ldsc.UKB_telomere_gwas_summarystats.N472174 \
	--a1 effec_allele \
	--a2 other_allele \
	--p p_value \
	--N 472174 \
	--snp variant_id \
	--frq effect_allele_frequency \
	--merge-alleles /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/w_hm3.snplist 

## Prep MPN: 
## zcat /medpop/esp/mesbah/GWAS_CHIP/MPN_summary/MPN_metaGWAS_sumstats.tsv.gz | cut -f2-10 > /broad/hptmp/mesbah/ukb_chip/mpn-gwas/MPN_metaGWAS_sumstats.c9.tsv &
/medpop/esp2/mesbah/tools/ldsc/munge_sumstats.py \
	--sumstats /broad/hptmp/mesbah/ukb_chip/mpn-gwas/MPN_metaGWAS_sumstats.c9.tsv \
	--out /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/ldsc.MPN_metaGWAS_sumstats.ncas3797_ncontrl1152977_N1156774 \
	--a1 ALT \
	--a2 REF \
	--p pvalue \
	--N 1156774  \
	--snp RSID \
	--frq MAF \
	--merge-alleles /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/w_hm3.snplist 

### mCA GWAS Maryam et al. 
#chrom	pos	rsid	ref	alt	neg_log_pvalue	beta	stderr_beta	alt_allele_freq
# 1	14933	rs199856693	G	A	0.23457734190617724	0.028401	0.051687	0.019737
## In the UKB cohort, of 444,199 unrelated individuals without a known history of hematologic malignancy, 66,011 (14.9%) carried an mCA (15,350 autosomal) and 12,398 (3.2%) carried an expanded mCA clone, defined as an mCA mutation present in at least 10% of peripheral leukocytes (2,985 autosomal)
## sample prev: 12398/444199 = 0.02791091 ~0.03
## pop prev (ukb + mgb + BBJ + FINGEN): 12398/444199 + 1026/22461 + 1676/125541 + 9558/175690 = 0.03533574 (CUB #+ 177/871)
/medpop/esp2/mesbah/tools/ldsc/munge_sumstats.py \
	--sumstats /medpop/esp/mesbah/GWAS_CHIP/mCA/Maryam.expanded_mCA_summary_stats.N444199.gz \
	--out /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/ldsc.Maryam.mCA_summary_stats.N444199.samP03_popP035 \
	--a1 alt \
	--a2 ref \
	--p neg_log_pvalue \
	--N 444199 \
	--snp rsid \
	--frq alt_allele_freq \
	--merge-alleles /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/w_hm3.snplist

### CHIP data
# qsub -N getsum -wd /broad/hptmp/mesbah/ukb_chip/tmpdir /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/prep_MetaGWAS.summary_fromVCF.sh /broad/hptmp/mesbah/ukb_chip/sumGWAS 

## /medpop/esp2/mesbah/tools/ldsc/munge_sumstats.py --sumstats /broad/hptmp/mesbah/ukb_chip/sumGWAS/lifted_hg37.GWAMA.meta_CHIP.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.tsv.gz --out /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/ldsc.meta_CHIP.allsamples.eaf001_min2Studies --a1 EffectAllele --a2 OtherAllele --p Pvalue --N-col N --snp RSID --frq EAF

## wget -P /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/ https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
# bunzip2 /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/w_hm3.snplist.bz2

while read files 
do 
	/medpop/esp2/mesbah/tools/ldsc/munge_sumstats.py \
		--sumstats ${files} \
		--out /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/ldsc.$(basename ${files} ".tsv.gz") \
		--a1 EffectAllele \
		--a2 OtherAllele \
		--p Pvalue \
		--N-col N \
		--snp RSID \
		--frq EAF \
		--merge-alleles /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/w_hm3.snplist 

done < <(ls -lh /broad/hptmp/mesbah/ukb_chip/sumGWAS/lifted_hg37.*.eaf001_min2Studies.tsv.gz | awk '{print $NF}')


## /medpop/esp2/mesbah/tools/ldsc/ldsc.py --rg /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/ldsc.lifted_hg37.GWAMA.EUR_only_CHIP_GWAS.TopMedUKbbMGbbBioVU.eaf001_min2Studies.sumstats.gz,/broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/ldsc.UKB_telomere_gwas_summarystats.N472174.sumstats.gz --ref-ld-chr /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/eur_w_ld_chr/ --w-ld-chr /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/eur_w_ld_chr/ --samp-prev 0.05,nan --pop-prev 0.05,nan --out /broad/hptmp/mesbah/ukb_chip/h2/genCor/EUR_only_CHIP_GWAS.UKB_telomere.prev5

## /medpop/esp2/mesbah/tools/ldsc/ldsc.py --rg /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/ldsc.lifted_hg37.GWAMA.EUR_only_CHIP_GWAS.TopMedUKbbMGbbBioVU.eaf001_min2Studies.sumstats.gz,/broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/ldsc.MPN_metaGWAS_sumstats.ncas3797_ncontrl1152977_N1156774.sumstats.gz --ref-ld-chr /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/eur_w_ld_chr/ --w-ld-chr /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/eur_w_ld_chr/ --samp-prev 0.05,0.0033 --pop-prev 0.05,2.17e-05 --out /broad/hptmp/mesbah/ukb_chip/h2/genCor/EUR_only_CHIP_GWAS.MPN_metaGWAS.prev5

