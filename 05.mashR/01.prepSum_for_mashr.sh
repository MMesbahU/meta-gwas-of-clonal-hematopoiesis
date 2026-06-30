#!/bin/bash/

### Prepare GWAS summary for mashR
mkdir -p /broad/hptmp/mesbah/gwas/mashr

#### MPN GWAS (Bao, E. L. et al. 2020)
# zcat /medpop/esp/mesbah/GWAS_CHIP/MPN_summary/MPN_metaGWAS_sumstats.tsv.gz | head -2
# MarkerName	RSID	CHR	POS	REF	ALT	Effect	StdErr	pvalue	MAF
# 1:768448_G_A	rs12562034	1	768448	G	A	-0.0253	0.0615	0.681	0.1076
## # Effect allele = ALT allele in the summary
## frequency of the minor allele is provided (instead of effect allele) 
zcat /medpop/esp/mesbah/GWAS_CHIP/MPN_summary/MPN_metaGWAS_sumstats.tsv.gz | cut -f1,7-10 | sed 's:_:\::g' | awk '(NR==1){print "SNPID\tBETA\tSE\tP\tMAF\tMPN_Z\tMPN_N"}(NR>1 && $1~/^[0-9]+/){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$2/$3"\t1156774"}' | awk '!seen[$1]++' | gzip -c > /broad/hptmp/mesbah/gwas/mashr/chr1_22.mpn_hg37_summary.tsv.gz &

#### LTL GWAS in UKB: Codd, V. et al. Polygenic basis and biomedical consequences of telomere length variation. Nat Genet 53, 1425-1433 (2021).  
# zcat /medpop/esp/mesbah/GWAS_CHIP/ltl_gwas/UKB_telomere_gwas_summarystats.tsv.gz | head -2                                                       
# variant_id	p_value	chromosome	base_pair_location	effec_allele	other_allele	effect_allele_frequency	beta	standard_error
# rs367896724	0.13	1	10177	A	AC	0.599214	0.00452419	0.00296934
## REF allele is the Effect allele
## converted to ALT allele by multiplying -1 * BETA
## alse changed Effect allele by (1-EAF)
zcat /medpop/esp/mesbah/GWAS_CHIP/ltl_gwas/UKB_telomere_gwas_summarystats.tsv.gz | awk '(NR==1){print "SNPID\tBETA\tSE\tP\tEAF\tLTL_Z\tLTL_N"}(NR>1 && $3~/^[0-9]+/){print $3":"$4":"$5":"$6"\t"(-1 * $8)"\t"$9"\t"$2"\t"(1-$7)"\t"((-1*$8) /$9)"\t472174"}'| awk '!seen[$1]++' | gzip -c > /broad/hptmp/mesbah/gwas/mashr/chr1_22.ltl_hg37_summary.tsv.gz &

##### Expanded mCA GWAS in UKB 500k
### mCA GWAS in UKB: Zekavat, S.M. et al. Hematopoietic mosaic chromosomal alterations increase the risk for diverse types of infection. Nat Med 27, 1012-1024 (2021)
# zcat /medpop/esp/mesbah/GWAS_CHIP/mCA/Maryam.expanded_mCA_summary_stats.N444199.gz | head -2
# chrom	pos	rsid	ref	alt	neg_log_pvalue	beta	stderr_beta	alt_allele_freq
# 1	14933	rs199856693	G	A	0.23457734190617724	0.028401	0.051687	0.019737
zcat /medpop/esp/mesbah/GWAS_CHIP/mCA/Maryam.expanded_mCA_summary_stats.N444199.gz | awk '(NR==1){print "SNPID\tBETA\tSE\tP\tEAF\tmCA_Z\tmCA_N"}(NR>1 && $1~/^[0-9]+/){print $1":"$2":"$4":"$5"\t"$7"\t"$8"\t"10^-($6)"\t"$9"\t"($7/$8)"\t444199"}' | awk '!seen[$1]++' | gzip -c > /broad/hptmp/mesbah/gwas/mashr/chr1_22.expanded_mCA_hg37_summary.tsv.gz & 


## Overall CHIP GWAS in 650k multi-ancestry samples
# zcat lifted_hg37.eur_metaGWAS.CHIP.GWAMA.hg37.eaf001_min2Studies.tsv.gz | head -2
# CHR	POS	hg38MarkerID	REF	ATL	EffectAllele	OtherAllele	EAF	BETA	SE	P	NDirection	Het_P
# 1	54490	chr1:54490:G:A	G	A	A	G	0.165455	-0.011681	0.028213	0.678886	406991	?-?+?	0.306857
while read trait;

do
	zcat /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/multi_ancestry_summary/lifted_hg37.metaGWAS.${trait}.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.hg37.eaf001_min2Studies.tsv.gz | awk -v mytrait=${trait} '(NR==1){print "SNPID\tBETA\tSE\tP\tEAF\t"mytrait"_Z\t"mytrait"_N"}(NR>1 && $1~/^[0-9]+/){print $1":"$2":"$4":"$5"\t"$9"\t"$10"\t"$11"\t"$8"\t"($9/$10)"\t"$12}' | awk '!seen[$1]++' | gzip -c > /broad/hptmp/mesbah/gwas/mashr/chr1_22.${trait}_hg37_summary.tsv.gz

done < <(echo -e "CHIP\nDNMT3A\nTET2") &

## DNMT3A CHIP
# zcat lifted_hg37.metaGWAS.DNMT3A.GWAMA.hg37.eaf001_min2Studies.tsv.gz | awk '(NR==1){print "SNPID\tBETA\tSE\tP\tEAF\tDNMT3A_Z\tDNMT3A_N"}(NR>1 && $1~/^[0-9]+/){print $1":"$2":"$4":"$5"\t"$9"\t"$10"\t"$11"\t"$8"\t"($9/$10)"\t"$12}' | awk '!seen[$1]++' | gzip -c > /broad/hptmp/mesbah/gwas/mashr/chr1_22.DNMT3A_hg37_summary.tsv.gz &

## TET2 CHIP
# zcat lifted_hg37.metaGWAS.TET2.GWAMA.hg37.eaf001_min2Studies.tsv.gz | awk '(NR==1){print "SNPID\tBETA\tSE\tP\tEAF\tTET2_Z\tTET2_N"}(NR>1 && $1~/^[0-9]+/){print $1":"$2":"$4":"$5"\t"$9"\t"$10"\t"$11"\t"$8"\t"($9/$10)"\t"$12}' | awk '!seen[$1]++' | gzip -c > /broad/hptmp/mesbah/gwas/mashr/chr1_22.TET2_hg37_summary.tsv.gz &

