#!/bin/bash

######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
#########
## Aug 8, 2022
# while read traits; do qsub -R y -wd /broad/hptmp/mesbah/gwas/ukb450k/ukb250vs200k -pe smp 1 -binding linear:1 -l h_vmem=30G -l h_rt=10:00:00 -N metaGWAS.${traits}.ukb200k_vs_250k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/01.run.GWAMA.ukb200k_vs_ukb250k.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/01.PrepSummary/gwama_input.has${traits}.UKB200k_vs_UKB250k.hg37.txt /broad/hptmp/mesbah/gwas/ukb450k/ukb250vs200k/metaGWAS.${traits}.ukb200k_vs_ukb250k; done < <(echo -e "CHIP\nDNMT3A\nTET2\nASXL1")
##########################
###
gwas_list=${1}

output_prefix=${2}

## Run GWAMA
## Name	Chr	Pos_hg19	Ref	Alt	Trait	Cohort	Model	Effect	LCI_Effect	UCI_Effect	Pval	AAF	Num_Cases	Cases_Ref	Cases_Het	Cases_Alt	Num_Controls	Controls_Ref	Controls_Het	Controls_Alt	REGENIE_BETA	REGENIE_SE	INFO	MAC	N
GWAMA=/medpop/esp2/mesbah/tools/GWAMA_v2.2.2/GWAMA 
${GWAMA} \
	-i ${gwas_list} \
	-qt \
	--name_marker Name \
	--name_n N \
	--name_ea Alt \
	--name_nea Ref \
	--name_eaf AAF \
	--name_beta REGENIE_BETA \
	--name_se REGENIE_SE \
	--indel_alleles \
	-o ${output_prefix}

## compress
gzip ${output_prefix}.out

### Extart SNP present in >=2 studies, MAF>=0.1%
# zcat  metaGWAS.CHIP.TOPMed2019UKB200kMGBB53kBioVU54k.out.gz | awk '(NR==1){print "RSID\tREF\tALT\tAAF\tBETA\tSE\tP"}(NR>1 && $1~/^chr[0-9]+/ && $15>1 && $4>=0.001 && $4<=0.999){print $1"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6"\t"$10}' | bgzip -c > chr1_22.n2maf001.CHIP.tsv.gz &
# zcat  metaGWAS.DNMT3A.TOPMed2019UKB200kMGBB53kBioVU54k.out.gz | awk '(NR==1){print "RSID\tREF\tALT\tAAF\tBETA\tSE\tP"}(NR>1 && $1~/^chr[0-9]+/ && $15>1 && $4>=0.001 && $4<=0.999){print $1"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6"\t"$10}' | bgzip -c > chr1_22.n2maf001.DNMT3A.tsv.gz &
# zcat  metaGWAS.TET2.TOPMed2019UKB200kMGBB53kBioVU54k.out.gz | awk '(NR==1){print "RSID\tREF\tALT\tAAF\tBETA\tSE\tP"}(NR>1 && $1~/^chr[0-9]+/ && $15>1 && $4>=0.001 && $4<=0.999){print $1"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6"\t"$10}' | bgzip -c > chr1_22.n2maf001.TET2.tsv.gz &
# zcat chr1_22.n2maf001.CHIP.tsv.gz | head -1 > sorted.chr1_22.n2maf001.CHIP.tsv && zcat chr1_22.n2maf001.CHIP.tsv.gz | awk 'NR>1' | sort -k1 -V >> sorted.chr1_22.n2maf001.CHIP.tsv && bgzip sorted.chr1_22.n2maf001.CHIP.tsv &
###########


## convert to vcf
# rs_number	reference_allele	other_allele	eaf	beta	se	beta_95L	beta_95U	z	p-value	_-log10_p-value	q_statistic	q_p-value	i2	n_studies	n_samples	effects
# chr1:13289:CCT:C	C	CCT	0.002043	0.310031	0.348919	-0.373849	0.993912	0.888549	0.374242	0.426848	0.000000	1.000000	-nan	1	29373	+???
#echo -e '##fileformat=VCFv4.2' > ${tmpVCF}
#echo -e '##INFO=<ID=rs_number,Number=1,Type=String,Description="hg38 chr:pos:ref:alt ids">'
#echo -e '##INFO=<ID=reference_allele,Number=1,Type=String,Description="GWAMA effect allele (reference allele for effect)">'
#echo -e '##INFO=<ID=other_allele,Number=1,Type=String,Description="GWAMA other allele (non reference allele for effect)">'
#echo -e '##INFO=<ID=eaf,Number=A,Type=Float,Description="GWAMA effect allele (reference allele for effect)">'
#echo -e '##INFO=<ID=reference_allele,Number=1,Type=String,Description="GWAMA effect allele (reference allele for effect)">'
#echo -e '##INFO=<ID=reference_allele,Number=1,Type=String,Description="GWAMA effect allele (reference allele for effect)">'

###INFO=<ID=AAF,Number=A,Type=Float,Description="Alternate Allele Frequency">\n##INFO=<ID=NCASES,Number=1,Type=Integer,Description="Number of CHIP Cases">\n##INFO=<ID=NCONTROL,Number=1,Type=Integer,Description="Number of CHIP Controls">\n##INFO=<ID=REGENIE_BETA,Number=A,Type=Float,Description="SNP effect estimate beta">\n##INFO=<ID=REGENIE_SE,Number=A,Type=Float,Description="SE estimate of SNP effect">\n##INFO=<ID=INFO,Number=A,Type=Float,Description="Imputation r square">\n##INFO=<ID=MAC,Number=1,Type=Integer,Description="Minor allele count">\n##INFO=<ID=N,Number=A,Type=Float,Description="Total Samples with genotype">\n##INFO=<ID=old_hg19_id,Number=1,Type=String,Description="UKB HG19/hg37 coordinates">\n##INFO=<ID=old_ukbrsid,Number=1,Type=String,Description="UKB rsids">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' > ${tmpVCF}

#zcat ${inGWAS} | awk 'NR>1{print $2"\t"$3"\t"$1"\t"$4"\t"$5"\t.\t.\tPVAL="$12";AAF="$13";REGENIE_BETA="$22";REGENIE_SE="$23";INFO="$24";MAC="$25";NCASES="$14";NCONTROL="$18";old_hg19_id="$2":"$3":"$4":"$5";old_ukbrsid="$1";N="($14+$18)}' >> ${tmpVCF}

#### 

######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

