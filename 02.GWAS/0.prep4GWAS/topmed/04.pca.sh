#!/bin/bash

## qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/tmpdir -N PC_APPROX -l h_vmem=20G -l h_rt=5:00:00 -pe smp 2 -binding linear:2 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/0.prep4GWAS/topmed/04.pca.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/plink/topmed_N74k.MultiANC.freeze.10b.chr1_22.plink_maf01hwe1e50geno10mind10 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/plink/MultiANC/qcdGSAsnp.topmed.MultiANC.snplist /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/pheno/qcd_topmed.PC20apx /medpop/esp2/mesbah/tools/plink2_linux_x86_64_20200831

plink_file=${1}

SNP_list=${2}

outFilePrefix=${3} # top20_approx_PC_qc_arrays_SNPs

PLINK2=${4} 

### Get top 20 PCs
### use Approximate approach
## exclude indels for PCA
## PCA Analysis
${PLINK2} \
	--bfile ${plink_file} \
	--extract <( sed 's:\::\t:g' ${SNP_list} | awk '{snp=$1":"$2":"$3":"$4;REF=$3;ALT=$4}( (REF=="A"||REF=="C"||REF=="G"||REF=="T") && (ALT=="A"||ALT=="C"||ALT=="G"||ALT=="T")){print snp}' ) \
	--pca 20 approx \
	--maf 0.05 \
	--hwe 1e-20 \
	--threads $(nproc --all) \
	--out ${outFilePrefix} \
	--memory 120000

