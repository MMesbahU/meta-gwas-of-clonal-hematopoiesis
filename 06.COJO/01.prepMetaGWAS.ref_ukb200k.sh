#!/bin/bash

# qsub -wd /broad/hptmp/mesbah/gwas/cojo/tmpdir -R y -l h_vmem=20G -l h_rt=20:00:00 -pe smp 1 -binding linear:1 -N sum_cojo /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/06.COJO/01.prepMetaGWAS.ref_ukb200k.sh summary_4_cojo /broad/hptmp/mesbah/gwas/cojo/input

# qsub -wd /broad/hptmp/mesbah/ukb_chip/tmpdir -R y -l h_vmem=20G -l h_rt=20:00:00 -pe smp 1 -binding linear:1 -N COJO_LDAK_prep /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.COJO/01.prepMetaGWAS.ref_ukb200k.sh in_cojo_eaf001_min2Studies.8March2022 in_ldak_eaf001_min2Studies.8March2022 /broad/hptmp/mesbah/ukb_chip/cojo/cojo_input /broad/hptmp/mesbah/ukb_chip/h2/input_ldak

source /broad/software/scripts/useuse

use Bcftools

# 
# MAF 0.1%
# n Studies>=2
# 
#
cojo_input=${1}

#ldak_input=${2}

cojoINdir=${2}

#ldakINdir=${4}

mkdir -p ${cojoINdir}

#mkdir -p ${ldakINdir}

for files in $(ls -lhv /medpop/esp/mesbah/GWAMA_VCF/hg37/lifted_hg37.*.hg37_dbSNP.vcf.gz | awk '{print $NF}')
do
	gwas_summary_v1=${cojoINdir}/ref_ukb200k.${cojo_input}.$(basename ${files} ".vcf.gz").tsv

	# gwas_summary_v2=${cojoINdir}/ref_mgbb53k.${cojo_input}.$(basename ${files} ".vcf.gz").tsv

	# Fileldak_input=${ldakINdir}/${ldak_input}.$(basename ${files} ".vcf.gz").txt

#### COJO input from VCF file
# COJO input fields: 
	echo -e "SNP\tA1\tA2\tfreq\tb\tse\tp\tN" > ${gwas_summary_v1}

	#
	bcftools view  -i 'EAF>=0.001 & EAF<=0.999 & N_Studies>=2' ${files} | bcftools query -f '%ID\t%INFO/EffectAllele\t%INFO/OtherAllele\t%INFO/EAF\t%INFO/BETA\t%INFO/SE\t%INFO/Pvalue\t%INFO/N_METAL\n' | awk '!seen[$1]++' >> ${gwas_summary_v1}
	
	## MGB Ref
	# echo -e "SNP\tA1\tA2\tfreq\tb\tse\tp\tN" > ${gwas_summary_v2}
	
	# bcftools view  -i 'EAF>=0.001 & EAF<=0.999 & N_Studies>=2' ${files} | bcftools query -f '%INFO/varID_hg38\t%INFO/EffectAllele\t%INFO/OtherAllele\t%INFO/EAF\t%INFO/BETA\t%INFO/SE\t%INFO/Pvalue\t%INFO/N_METAL\n' | awk '!seen[$1]++' >> ${gwas_summary_v2}
# LDAK input fields: "Predictor A1 A2 Z n"
# only High quality SNP
## There is Z stats also from GWAMA: %INFO/Z
#	bcftools view  -i 'EAF>=0.001 & EAF<=0.999 & N_Studies>=2' ${files} | bcftools query -f '%CHROM %POS %INFO/EffectAllele %INFO/other_allele %INFO/BETA %INFO/SE %INFO/N_METAL\n' | awk '{snp=$1":"$2;a1=$3;a2=$4;z=($5/$6);n=$7}(NR==1){print "Predictor A1 A2 Z n"}(NR>1 && (a1=="A"||a1=="C"||a1=="G"||a1=="T") && (a2=="A"||a2=="C"||a2=="G"||a2=="T")){print snp, a1, a2, z, n}' | awk '!seen[$1]++' > ${Fileldak_input}

done 


