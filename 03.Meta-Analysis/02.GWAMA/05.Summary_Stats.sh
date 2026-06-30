#!/bin/bash

##### Extract Variants from VCF files
## EAF>=1% & P<5e-3
# qsub -R y -wd /broad/hptmp/mesbah/gwas/lifted_hg37 -N meta_summary -pe smp 1 -binding linear:1 -l h_rt=20:00:00 -l h_vmem=20G /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/05.Summary_Stats.sh

source /broad/software/scripts/useuse
use Bcftools
use Tabix
### Nov 15, 2021
# N_Studies >=2
# MAF>=0.1%

### Pvalue contains '0' in many instances; have to use 10^-log10_Pvalue
## convert log10(pvalue) to pvalue

# for files in $(ls -lhrt /medpop/esp/mesbah/GWAMA_VCF/hg37/lifted_hg37.eur_metaGWAS.*.GWAMA.hg37_dbSNP.vcf.gz | awk '{print $NF}'); do echo -e "CHR\tPOS\tMarkerID\tREF\tATL\tEffectAllele\tOtherAllele\tEAF\tBETA\tSE\tP\tN\tDirection\tHet_P" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/eur_only_summary/$(basename ${files} ".vcf.gz").eaf001_min2Studies.tsv && bcftools view  -i 'EAF>=0.001 & EAF<=0.999 & N_Studies>=2' ${files} | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/EffectAllele\t%INFO/OtherAllele\t%INFO/EAF\t%INFO/BETA\t%INFO/SE\t%INFO/log10_Pvalue\t%INFO/N_METAL\t%INFO/Effect_Direction\t%INFO/Q_Pvalue\n' | awk '$1~/^[0-9]+/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"10^-$11"\t"$12"\t"$13"\t"$14}' >> /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/eur_only_summary/$(basename ${files} ".vcf.gz").eaf001_min2Studies.tsv && bgzip -f /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/eur_only_summary/$(basename ${files} ".vcf.gz").eaf001_min2Studies.tsv; done &
##
for files in $(ls -lhrt /medpop/esp/mesbah/GWAMA_VCF/hg37/lifted_hg37.eur_metaGWAS.*.GWAMA.hg37_dbSNP.vcf.gz | awk '{print $NF}')

do 
	# rsid or hg38MarkerID in Column 3
	echo -e "CHR\tPOS\tMarkerID\tREF\tATL\tEffectAllele\tOtherAllele\tEAF\tBETA\tSE\tP\tN\tDirection\tHet_P" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/eur_only_summary/$(basename ${files} ".vcf.gz").eaf001_min2Studies.tsv

#
	bcftools view  -i 'EAF>=0.001 & EAF<=0.999 & N_Studies>=2' ${files} | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/EffectAllele\t%INFO/OtherAllele\t%INFO/EAF\t%INFO/BETA\t%INFO/SE\t%INFO/log10_Pvalue\t%INFO/N_METAL\t%INFO/Effect_Direction\t%INFO/Q_Pvalue\n' | awk '$1~/^[0-9]+/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"10^-$11"\t"$12"\t"$13"\t"$14}' >> /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/eur_only_summary/$(basename ${files} ".vcf.gz").eaf001_min2Studies.tsv

	bgzip -f /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/eur_only_summary/$(basename ${files} ".vcf.gz").eaf001_min2Studies.tsv
	# need sorted file for index
	# tabix -s 1 -b 2 -e 2 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/eur_only_summary/$(basename ${files} ".vcf.gz").eaf001_min2Studies.tsv.gz

done

## Multi-Ancestry GWAS
## for files in $(ls -lhrt /medpop/esp/mesbah/GWAMA_VCF/hg37/lifted_hg37.metaGWAS.*.hg37_dbSNP.vcf.gz | awk '{print $NF}'); do echo -e "CHR\tPOS\tMarkerID\tREF\tALT\tEffectAllele\tOtherAllele\tEAF\tBETA\tSE\tP\tN\tDirection\tHet_P" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/multi_ancestry_summary/$(basename ${files} ".vcf.gz").eaf001_min2Studies.tsv && bcftools view  -i 'EAF>=0.001 & EAF<=0.999 & N_Studies>=2' ${files} | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/EffectAllele\t%INFO/OtherAllele\t%INFO/EAF\t%INFO/BETA\t%INFO/SE\t%INFO/log10_Pvalue\t%INFO/N_METAL\t%INFO/Effect_Direction\t%INFO/Q_Pvalue\n' | awk '$1~/^[0-9]+/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"10^-$11"\t"$12"\t"$13"\t"$14}' >> /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/multi_ancestry_summary/$(basename ${files} ".vcf.gz").eaf001_min2Studies.tsv && bgzip -f /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/multi_ancestry_summary/$(basename ${files} ".vcf.gz").eaf001_min2Studies.tsv; done &
# /medpop/esp/mesbah/GWAS_CHIP/meta_gwas/meta_uk450k/450k/hg38/metaGWAS.*hg38.vcf.gz
##
for files in $(ls -lhrt /medpop/esp/mesbah/GWAMA_VCF/hg37/lifted_hg37.metaGWAS.*.hg37_dbSNP.vcf.gz | awk '{print $NF}')

do
	# rsid or hg38MarkerID in Column 3
        echo -e "CHR\tPOS\tMarkerID\tREF\tALT\tEffectAllele\tOtherAllele\tEAF\tBETA\tSE\tP\tN\tDirection\tHet_P" > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/multi_ancestry_summary/$(basename ${files} ".vcf.gz").eaf001_min2Studies.tsv

##
	bcftools view  -i 'EAF>=0.001 & EAF<=0.999 & N_Studies>=2' ${files} | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/EffectAllele\t%INFO/OtherAllele\t%INFO/EAF\t%INFO/BETA\t%INFO/SE\t%INFO/log10_Pvalue\t%INFO/N_METAL\t%INFO/Effect_Direction\t%INFO/Q_Pvalue\n' | awk '$1~/^[0-9]+/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"10^-$11"\t"$12"\t"$13"\t"$14}' >> /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/multi_ancestry_summary/$(basename ${files} ".vcf.gz").eaf001_min2Studies.tsv

	bgzip -f /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/multi_ancestry_summary/$(basename ${files} ".vcf.gz").eaf001_min2Studies.tsv

	# tabix -s 1 -b 2 -e 2 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/multi_ancestry_summary/$(basename ${files} ".vcf.gz").eaf001_min2Studies.tsv.gz

done

#### P < 5e-8
# while read trait; do echo -e "CHR\tPOS\tMarkerID\tvarID_hg38\tREF\tATL\tEffectAllele\tOtherAllele\tEAF\tBETA\tSE\tP\tN\tDirection\tHet_P" > multiancestry.${trait}.p5e8_eaf001_min2Studies.tsv && bcftools view  -i 'EAF>=0.001 & EAF<=0.999 & N_Studies>=2 & Pvalue<5e-8' /medpop/esp/mesbah/GWAMA_VCF/hg38/metaGWAS.${trait}.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.hg38_dbSNP.vcf.gz | bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/varID_hg38\t%REF\t%ALT\t%INFO/EffectAllele\t%INFO/OtherAllele\t%INFO/EAF\t%INFO/BETA\t%INFO/SE\t%INFO/Pvalue\t%INFO/N_METAL\t%INFO/Effect_Direction\t%INFO/Q_Pvalue\n'  >> multiancestry.${trait}.p5e8_eaf001_min2Studies.tsv; done < <(echo -e "CHIP\nDNMT3A\nTET2") &

# while read trait; do echo -e "CHR\tPOS\tMarkerID\tvarID_hg38\tREF\tATL\tEffectAllele\tOtherAllele\tEAF\tBETA\tSE\tP\tN\tDirection\tHet_P" > eur_only.${trait}.p5e8_eaf001_min2Studies.tsv && bcftools view  -i 'EAF>=0.001 & EAF<=0.999 & N_Studies>=2 & Pvalue<5e-8' /medpop/esp/mesbah/GWAMA_VCF/hg38/eur_metaGWAS.${trait}.GWAMA.hg38_dbSNP.vcf.gz | bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/varID_hg38\t%REF\t%ALT\t%INFO/EffectAllele\t%INFO/OtherAllele\t%INFO/EAF\t%INFO/BETA\t%INFO/SE\t%INFO/Pvalue\t%INFO/N_METAL\t%INFO/Effect_Direction\t%INFO/Q_Pvalue\n' >> eur_only.${trait}.p5e8_eaf001_min2Studies.tsv; done < <(echo -e "CHIP\nDNMT3A\nTET2") &



