#!/bin/bash

## qsub -wd /broad/hptmp/mesbah/gwas/sumher -R y -l h_vmem=20G -l h_rt=20:00:00 -pe smp 1 -binding linear:1 -N ldak_input /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.SNP_heritability/01.prepMetaGWAS.sh /broad/hptmp/mesbah/gwas/sumher /broad/hptmp/mesbah/gwas/sumher

## qsub -wd /broad/hptmp/mesbah/ukb_chip/tmpdir -R y -l h_vmem=20G -l h_rt=20:00:00 -pe smp 1 -binding linear:1 -N LDAK_prep /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/09.Heritability/01.prepMetaGWAS.sh in_cojo_eaf001_min2Studies.8March2022 in_ldak_eaf001_min2Studies.8March2022 /broad/hptmp/mesbah/ukb_chip/cojo/cojo_input /broad/hptmp/mesbah/ukb_chip/h2/input_ldak

source /broad/software/scripts/useuse

use Bcftools

# 
# MAF 0.1%
# n Studies>=2
# 
#
# cojo_input=${1}

ldak_input=${1}

# cojoINdir=${3}

ldakINdir=${2}

# mkdir -p ${cojoINdir}

mkdir -p ${ldakINdir}

for files in $(ls -lhv /medpop/esp/mesbah/GWAMA_VCF/hg37/lifted_hg37.*.hg37_dbSNP.vcf.gz | awk '{print $NF}')
do
#	Filecojo_input=${cojoINdir}/${cojo_input}.$(basename ${files} ".vcf.gz").tsv

	Fileldak_input=${ldakINdir}/${ldak_input}.$(basename ${files} ".vcf.gz").txt

#### COJO input from VCF file
# COJO input fields: echo -e "SNP\tA1\tA2\tfreq\tb\tse\tp\tN" > ${cojo_input} 
#	bcftools view  -i 'EAF>=0.001 & EAF<=0.999 & N_Studies>=2' ${files} | bcftools query -f '%ID\t%INFO/EffectAllele\t%INFO/other_allele\t%INFO/EAF\t%INFO/BETA\t%INFO/SE\t%INFO/log10_pvalue\t%INFO/N_METAL\n' | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"10^-$7"\t"$8}' | awk '{SNP=$1;A1=$2;A2=$3;freq=$4;b=$5;se=$6;p=$7;N=$8}(NR==1){print "SNP\tA1\tA2\tfreq\tb\tse\tp\tN"}(NR>1){print SNP"\t"A1"\t"A2"\t"freq"\t"b"\t"se"\t"p"\t"N}' | awk '!seen[$1]++' > ${Filecojo_input}

# LDAK input fields: "Predictor A1 A2 Z n"
# only High quality SNP
## There is Z stats also from GWAMA: %INFO/Z
	echo -e "Predictor A1 A2 Z n" > ${Fileldak_input}

	bcftools view  -i 'EAF>=0.001 & EAF<=0.999 & N_Studies>=2' ${files} | bcftools query -f '%CHROM %POS %INFO/EffectAllele %INFO/OtherAllele %INFO/Z %INFO/N_METAL\n' | awk '{snp=$1":"$2;a1=$3;a2=$4;z=$5;n=$6}($1~/^[0-9]+/ && (a1=="A"||a1=="C"||a1=="G"||a1=="T") && (a2=="A"||a2=="C"||a2=="G"||a2=="T")){print snp, a1, a2, z, n}' | awk '!seen[$1]++' >> ${Fileldak_input}

done 


