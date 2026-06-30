#!/bin/bash

#
## 2024 April:
## while read ANC; do qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/tmpdir -N qc_plink_${ANC} -l h_vmem=20G -l h_rt=5:00:00 -pe smp 2 -binding linear:2 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/0.prep4GWAS/topmed/03.qc_geno.sh /medpop/esp2/mesbah/tools/plink2_linux_x86_64_20200831 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/plink/topmed_N74k.MultiANC.freeze.10b.chr1_22.plink_maf01hwe1e50geno10mind10 1e-30 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/pheno/TOPMed_GWAS_pheno.${ANC}.04_23_2024.tsv.gz /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/plink qcdGSAsnp.topmed.${ANC} ${ANC}; done < <(echo -e "Female\nMale\nAFR\nAMR\nEUR\nMultiANC") 

## while read ANC; do qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/tmpdir -N qc_plink_${ANC} -l h_vmem=20G -l h_rt=5:00:00 -pe smp 2 -binding linear:2 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/0.prep4GWAS/topmed/03.qc_geno.sh /medpop/esp2/mesbah/tools/plink2_linux_x86_64_20200831 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/plink/topmed_N74k.MultiANC.freeze.10b.chr1_22.plink_maf01hwe1e50geno10mind10 1e-10 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/pheno/TOPMed_GWAS_pheno.${ANC}.04_23_2024.v02.tsv.gz /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/plink qcdGSAsnp.topmed.${ANC} ${ANC}; done < <(echo -e "EAS\nSAS")

PLINK2=${1}

PLINK_Prefix=${2}

HWE_P=${3}

Samples2Keep=${4}

outDir=${5}

OutPrefix=${6}

ANC=${7}

###
mkdir -p ${outDir}/${ANC}

### Plink filter
${PLINK2} \
	--bfile ${PLINK_Prefix} \
	--maf 0.01 \
	--geno 0.1 \
	--hwe ${HWE_P} \
	--mind 0.1 \
	--keep <(zcat ${Samples2Keep} | awk 'NR>1{print $1,$2}' ) \
	--write-snplist \
	--write-samples \
	--no-id-header \
	--threads $(nproc --all) \
	--out ${outDir}/${ANC}/${OutPrefix}



