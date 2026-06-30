#!/bin/bash

## for files in $(ls -l /medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/chr1_22.has*.21Aug2021_ukbEUR.tsv.gz | awk '{print $NF}'); do qsub -R y -wd /broad/hptmp/mesbah/ukb_chip/tmpdir -N prepSum4ldak -l h_rt=20:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/09.Heritability/0.run_prepSum4ldak.sh ${files} /broad/hptmp/mesbah/ukb_chip/h2/input_ldak; done
## wget https://genetics.ghpc.au.dk/doug/bld.ldak.hapmap.gbr.tagging.gz
## gzip -d bld.ldak.hapmap.gbr.tagging.gz
#### 
# source /broad/software/scripts/useuse

# use Bcftools
# use Tabix


### Read command line inputs
gwas_tsv=${1}

ldak_input_dir=${2}

### Prepare COJO input from VCF file
ldak_input=${ldak_input_dir}/summary_for_ldak.$(basename ${gwas_tsv} ".gz")

# echo -e "Predictor A1 A2 Z n" > ${gwas_summary}

## http://dougspeed.com/summary-statistics/
## info>=0.95, MAF>=0.001

## UKB data
zcat ${gwas_tsv} | awk '(NR>1){snp=$2":"$3;a1=$5;a2=$4;z=($22/$23);n=($14+$18)}(NR==1){print "Predictor A1 A2 Z n"}(NR>1 && (a1=="A"||a1=="C"||a1=="G"||a1=="T") && (a2=="A"||a2=="C"||a2=="G"||a2=="T") && $24>=0.95 && $13>=0.001 && $13<=0.999){print snp, a1, a2, z, n}' | awk '!seen[$1]++' > ${ldak_input} 


### TOPMed Eur data
# 1	13289	rs538791886	CCT	C	.	PASS	AC_Allele2=120;AF_Allele2=0.00204269;Allele1=CCT;Allele2=C;BETA=0.310031;MarkerID=chr1_13289_CCT_C;N=29373;SE=0.348919;Tstat=2.54658;p.value=0.374246;p.value.NA=0.374246

## LDAK
# files=/medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/eur_topmed/topmed.hg37.CHIP.EUR.results.vcf.gz

# echo -e "Predictor A1 A2 Z n" > /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/summary_for_ldak.$(basename ${files} ".vcf.gz").txt

# bcftools view  -i 'AF_Allele2>=0.001 & AF_Allele2<=0.999' ${files} | bcftools query -f '%CHROM %POS %INFO/Allele2 %INFO/Allele1 %INFO/BETA %INFO/SE %INFO/N\n' | awk '{snp=$1":"$2;a1=$3;a2=$4;z=($5/$6);n=$7}(NR==1){print "Predictor A1 A2 Z n"}(NR>1 && (a1=="A"||a1=="C"||a1=="G"||a1=="T") && (a2=="A"||a2=="C"||a2=="G"||a2=="T")){print snp, a1, a2, z, n}' | awk '!seen[$1]++' > /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/summary_for_ldak.$(basename ${files} ".vcf.gz").txt

## LDSC
# echo -e "rsid\ta1\ta2\tn\tp\tb" > /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/summary_for_ldsc.$(basename ${files} ".vcf.gz").tsv

# bcftools view  -i 'AF_Allele2>=0.001 & AF_Allele2<=0.999' ${files} | bcftools query -f '%ID\t%INFO/Allele2\t%INFO/Allele1\t%INFO/N\t%INFO/p.value\t%INFO/BETA\n' | awk '{snp=$1;a1=$2;a2=$3;n=$4;p=$5;b=$6}(NR==1){print "rsid\ta1\ta2\tn\tp\tb"}(NR>1 && (a1=="A"||a1=="C"||a1=="G"||a1=="T") && (a2=="A"||a2=="C"||a2=="G"||a2=="T")){print snp"\t"a1"\t"a2"\t"n"\t"p"\t"b}' | awk '$1!="."' | awk '!seen[$1]++' > /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/summary_for_ldsc.$(basename ${files} ".vcf.gz").tsv

## Eur meta: 4 studies
# zcat /medpop/esp2/mesbah/Meta_GWAS/Sep2021/FineMap_GWAMA_EAF1/lifted_hg37.GWAMA.EUR_only_CHIP_GWAS.TopMedUKbbMGbbBioVU.eaf001_min2Studies.tsv.gz | head -2
# CHR	POS	REF	ALT	RSID	hg38MarkerID	EffectAllele	OtherAllele	EAF	BETA	SE	Pvalue	N	Direction	Het_P
# 1	64931	G	A	rs62639104	chr1:64931:G:A	A	G	0.075119	-0.045242	0.037539	0.228104	197086	--??	0.939827

# zcat /medpop/esp2/mesbah/Meta_GWAS/Sep2021/FineMap_GWAMA_EAF1/lifted_hg37.GWAMA.EUR_only_CHIP_GWAS.TopMedUKbbMGbbBioVU.eaf001_min2Studies.tsv.gz | awk '(NR>1){snp=$1":"$2;a1=$7;a2=$8;z=($10/$11);n=$13}(NR==1){print "Predictor A1 A2 Z n"}(NR>1 && (a1=="A"||a1=="C"||a1=="G"||a1=="T") && (a2=="A"||a2=="C"||a2=="G"||a2=="T") && $9>=0.001 && $9<=0.999){print snp, a1, a2, z, n}' | awk '!seen[$1]++' > /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/summary_for_ldak.lifted_hg37.GWAMA.EUR_only_CHIP_GWAS.TopMedUKbbMGbbBioVU.eaf001_min2Studies.txt

######### MPN
## Prep MPN GWAS
# zcat /medpop/esp/mesbah/GWAS_CHIP/MPN_summary/MPN_metaGWAS_sumstats.tsv.gz | head -2
# MarkerName	RSID	CHR	POS	REF	ALT	Effect	StdErr	pvalue	MAF
# 1:768448_G_A	rs12562034	1	768448	G	A	-0.0253	0.0615	0.681	0.1076

zcat /medpop/esp/mesbah/GWAS_CHIP/MPN_summary/MPN_metaGWAS_sumstats.tsv.gz | awk '(NR>1){snp=$3":"$4;a1=$6;a2=$5;z=($7/$8)}(NR==1){print "Predictor A1 A2 Z n"}(NR>1 && (a1=="A"||a1=="C"||a1=="G"||a1=="T") && (a2=="A"||a2=="C"||a2=="G"||a2=="T") && $10<=0.999 && $10>=0.001){print snp, a1, a2, z, (3797 + 1152977)}' | awk '!seen[$1]++' > /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/input_MPN.GWAS.txt &

# Prepare MPN data downloded from GWAS catalog
# /medpop/esp/mesbah/GWAS_CHIP/MPN_summary/GCST90000032_GRCh37.tsv.gz | head -2
# variant_id	UKID	p_value	chromosome	base_pair_location	effect_allele	other_allele	effect_allele_frequency	odds_ratio	beta	standard_error
# rs12562034	1:768448_G_A	0.236	1	768448	A	G	0.116	1.11215547310647	0.1063	0.0897
# zcat /medpop/esp/mesbah/GWAS_CHIP/MPN_summary/GCST90000032_GRCh37.tsv.gz | awk '(NR>1){snp=$4":"$5;a1=$6;a2=$7;z=($10/$11)}(NR==1){print "Predictor A1 A2 Z n"}(NR>1 && (a1=="A"||a1=="C"||a1=="G"||a1=="T") && (a2=="A"||a2=="C"||a2=="G"||a2=="T") && $8<=0.999 && $8>=0.001){print snp, a1, a2, z, (3797 + 1152977)}' | awk '!seen[$1]++' > /broad/hptmp/mesbah/ukb_chip/h2/input_ldak/input_MPN_GWAS.GCST90000032_GRCh37.txt &

########################### LTL GWAS summary from UKB
## Paper: https://pubmed.ncbi.nlm.nih.gov/34611362/ 
## Polygenic basis and biomedical consequences of telomere length variation
# wget https://figshare.com/ndownloader/files/28414941?private_link=caa99dc0f76d62990195
# mv 28414941\?private_link\=caa99dc0f76d62990195 UKB_telomere_gwas_summarystats.tsv.gz
# variant_id    p_value chromosome      base_pair_location      effec_allele    other_allele    effect_allele_frequencybeta     standard_error
# rs367896724   0.13    1       10177   A       AC      0.599214        0.00452419      0.00296934
## Sample size: UKB 500k = 472174
zcat /medpop/esp/mesbah/GWAS_CHIP/ltl_gwas/UKB_telomere_gwas_summarystats.tsv.gz | awk '(NR>1){snp=$3":"$4;a1=$5;a2=$6;z=($8/$9)}(NR==1){print "Predictor A1 A2 Z n"}(NR>1 && (a1=="A"||a1=="C"||a1=="G"||a1=="T") && (a2=="A"||a2=="C"||a2=="G"||a2=="T") && $7<=0.999 && $7>=0.001){print snp, a1, a2, z, 472174}' | awk '!seen[$1]++' > /broad/hptmp/mesbah/ukb_chip/ltl_gwas/in_ldak.UKB_telomere_gwas_summarystats.NatGen_N472174.GRCh37.txt &


########### Beta, Z , and P value correlation for CHIP P<1.67e-8 (GWAS significant for 3 traits gwas 5e-8 / 3) 
	## hg37
	## variant_id    p_value chromosome      base_pair_location      effec_allele    other_allele    effect_allele_frequency beta     standard_error
	# rs367896724   0.13    1       10177   A       AC      0.599214        0.00452419      0.00296934

	## LTL Summary
	## chr:pos:noneffect:effect noneffect_a effect_a rsid eaf beta se z pval N
# zcat /medpop/esp/mesbah/GWAS_CHIP/ltl_gwas/UKB_telomere_gwas_summarystats.tsv.gz | awk '(NR>1){snp=$3":"$4":"$6":"$5;noneffect_a=$6;effect_a=$5;rsid=$1;eaf=$7;beta=$8;se=$9;z=($8/$9);pval=$2}(NR==1){print "Predictor noneffectAllele effectAllele RSID EAF BETA SE Z Pval N"}(NR>1 && $7<=0.999 && $7>=0.001){print snp, noneffect_a, effect_a, rsid, eaf, beta, se, z, pval, 472174}' | awk '!seen[$1]++' > /broad/hptmp/mesbah/ukb_chip/ltl_gwas/in_correlation.UKB_telomere_gwas_summarystats.NatGen_N472174.GRCh37.txt.gz


