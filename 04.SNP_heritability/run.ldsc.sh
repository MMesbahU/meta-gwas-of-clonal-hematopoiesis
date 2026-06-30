#!/bin/bash

## /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/09.Heritability/run.ldsc.sh /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/summary_for_ldsc.topmed.hg37.CHIP.EUR.results.tsv /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/eur_w_ld_chr /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/eur_w_ld_chr /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/Her_ldsc.topmed.hg37.CHIP.EUR.results

source /broad/software/scripts/useuse

use Anaconda

# cd /medpop/esp2/mesbah/tools/ldsc

# conda env create --file environment.yml
source activate ldsc
# source deactivate

# https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation
gwas_summary=${1} # /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/summary_for_ldsc.topmed.hg37.CHIP.EUR.results.tsv

ref_ld_dir={2} # /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/eur_w_ld_chr

w_ld_dir=${3} # /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/eur_w_ld_chr

outHer=${4} # 
	## get LD for EUR 
	# wget -P /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/ https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
	# tar -jxvf /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/eur_w_ld_chr.tar.bz2

## data formating
# /medpop/esp2/mesbah/tools/ldsc/munge_sumstats.py --out ldsc.topmed.hg37.CHIP.EUR.results --sumstats /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/summary_for_ldsc.topmed.hg37.CHIP.EUR.results.tsv 

	## Run h2 in EUR
/medpop/esp2/mesbah/tools/ldsc/ldsc.py \
	--h2 ${gwas_summary} \
	--ref-ld-chr ${ref_ld_dir}/ \
	--w-ld-chr ${w_ld_dir}/ \
	--out ${outHer}

## /medpop/esp2/mesbah/tools/ldsc/ldsc.py --h2 /broad/hptmp/mesbah/ukb_chip/h2/ldsc.topmed.hg37.CHIP.EUR.results.sumstats.gz --ref-ld-chr /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/eur_w_ld_chr/ --w-ld-chr /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/eur_w_ld_chr/ --out /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/Her_ldsc.topmed.hg37.CHIP.EUR.results
## liability scale
## /medpop/esp2/mesbah/tools/ldsc/ldsc.py --h2 /broad/hptmp/mesbah/ukb_chip/h2/ldsc.topmed.hg37.CHIP.EUR.results.sumstats.gz --ref-ld-chr /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/eur_w_ld_chr/ --w-ld-chr /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/eur_w_ld_chr/ --out /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/Her_ldsc.topmed.hg37.CHIP.EUR.results --samp-prev 0.05 --pop-prev 0.05


### Prepare UKB CHIP GWAS
# /medpop/esp2/mesbah/tools/ldsc/munge_sumstats.py --sumstats /medpop/esp2/mesbah/projects/Meta_GWAS/Sep2021/Zhi_allUKB/chr1_22.hasCHIP.11Aug2021_ukb200k.allUKB_MAF001_INFO_03.tsv.gz --out /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/ldsc.chr1_22.hasCHIP.11Aug2021_ukb200k.allUKB_MAF001_INFO_03 --a1 EffectAllele --a2 OtherAllele --p Pvalue --N-col N --snp RSID --frq EAF

## Prep LTL GWAS: 
# /medpop/esp2/mesbah/tools/ldsc/munge_sumstats.py --sumstats /medpop/esp/mesbah/GWAS_CHIP/ltl_gwas/UKB_telomere_gwas_summarystats.tsv.gz --out /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/ldsc.UKB_telomere_gwas_summarystats.N472174 --a1 effec_allele --a2 other_allele --p p_value --N 472174 --snp variant_id --frq effect_allele_frequency --merge-alleles /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/w_hm3.snplist 

## Prep MPN: 
## zcat /medpop/esp/mesbah/GWAS_CHIP/MPN_summary/MPN_metaGWAS_sumstats.tsv.gz | cut -f2-10 > /broad/hptmp/mesbah/ukb_chip/mpn-gwas/MPN_metaGWAS_sumstats.c9.tsv &
# /medpop/esp2/mesbah/tools/ldsc/munge_sumstats.py --sumstats /broad/hptmp/mesbah/ukb_chip/mpn-gwas/MPN_metaGWAS_sumstats.c9.tsv --out /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/ldsc.MPN_metaGWAS_sumstats.ncas3797_ncontrl1152977_N1156774 --a1 ALT --a2 REF --p pvalue --N 1156774  --snp RSID --frq MAF --merge-alleles /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/w_hm3.snplist 

### CHIP data
# qsub -N getsum -wd /broad/hptmp/mesbah/ukb_chip/tmpdir /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/prep_MetaGWAS.summary_fromVCF.sh /broad/hptmp/mesbah/ukb_chip/sumGWAS 

## /medpop/esp2/mesbah/tools/ldsc/munge_sumstats.py --sumstats /broad/hptmp/mesbah/ukb_chip/sumGWAS/lifted_hg37.GWAMA.meta_CHIP.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.tsv.gz --out /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/ldsc.meta_CHIP.allsamples.eaf001_min2Studies --a1 EffectAllele --a2 OtherAllele --p Pvalue --N-col N --snp RSID --frq EAF

## wget -P /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/ https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
# bunzip2 /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/w_hm3.snplist.bz2

# while read files; do /medpop/esp2/mesbah/tools/ldsc/munge_sumstats.py --sumstats ${files} --out /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/ldsc.$(basename ${files} ".tsv.gz") --a1 EffectAllele --a2 OtherAllele --p Pvalue --N-col N --snp RSID --frq EAF --merge-alleles /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/w_hm3.snplist; done < <(ls -lh /broad/hptmp/mesbah/ukb_chip/sumGWAS/lifted_hg37.*.eaf001_min2Studies.tsv.gz | awk '{print $NF}')

## /medpop/esp2/mesbah/tools/ldsc/ldsc.py --rg /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/ldsc.lifted_hg37.GWAMA.EUR_only_CHIP_GWAS.TopMedUKbbMGbbBioVU.eaf001_min2Studies.sumstats.gz,/broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/ldsc.UKB_telomere_gwas_summarystats.N472174.sumstats.gz --ref-ld-chr /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/eur_w_ld_chr/ --w-ld-chr /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/eur_w_ld_chr/ --samp-prev 0.05,nan --pop-prev 0.05,nan --out /broad/hptmp/mesbah/ukb_chip/h2/genCor/EUR_only_CHIP_GWAS.UKB_telomere.prev5

## /medpop/esp2/mesbah/tools/ldsc/ldsc.py --rg /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/ldsc.lifted_hg37.GWAMA.EUR_only_CHIP_GWAS.TopMedUKbbMGbbBioVU.eaf001_min2Studies.sumstats.gz,/broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/ldsc.MPN_metaGWAS_sumstats.ncas3797_ncontrl1152977_N1156774.sumstats.gz --ref-ld-chr /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/eur_w_ld_chr/ --w-ld-chr /broad/hptmp/mesbah/ukb_chip/h2/input_ldsc/eur_w_ld_chr/ --samp-prev 0.05,0.0033 --pop-prev 0.05,2.17e-05 --out /broad/hptmp/mesbah/ukb_chip/h2/genCor/EUR_only_CHIP_GWAS.MPN_metaGWAS.prev5

