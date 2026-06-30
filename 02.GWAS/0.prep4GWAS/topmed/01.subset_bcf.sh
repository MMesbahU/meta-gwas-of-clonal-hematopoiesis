#!/bin/bash

## qsub -t 1-22 -R y -N subset_bcf -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/tmpdir -l h_rt=20:00:00 -l h_vmem=20G -pe smp 6 -binding linear:6 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/0.prep4GWAS/topmed/01.subset_bcf.sh /medpop/esp2/projects/topmed/freeze10/phased /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/pheno/ANC_KNN.chip_input_01_filter_04_19_2024.tsv /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/sort.chr1_22.hg38.GSA_24v3_0_A2.tsv topmed_N74k.MultiANC /medpop/esp2/mesbah/tools/plink2_linux_x86_64_20200831 5 

source /broad/software/scripts/useuse

use Bcftools

use Tabix

##
genodir=${1} #/medpop/esp2/projects/topmed/freeze10/phased
 
outDir=${2} #/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno

samples2keep=${3}

chr_pos2keep=${4}

outPrefix=${5} # topmed_N74k.MultiANC

plink2=${6}

CPU=${7}

###
chr=${SGE_TASK_ID}

bcf_file="${genodir}/freeze.10b.chr${chr}.pass_only.phased.bcf"

bcf_outfile="${outDir}/${outPrefix}.freeze.10b.chr${chr}.pass_only_phased.bcftools.vcf.gz"

plink_outfile="${outDir}/${outPrefix}.freeze.10b.chr${chr}.plink_maf01hwe1e30geno10mind10"


### BCFTOOLS
# no header in snp file 
bcftools view ${bcf_file} \
	-S <(awk 'NR>1 {print $1}' ${samples2keep} ) \
	-R <( awk '{print $1"\t"$2}' ${chr_pos2keep} ) \
	-Oz \
	--threads ${CPU} \
	-o ${bcf_outfile}
#### add chr:pos:ref:alt
# for myfile in $(ls -lh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/topmed_N74k.MultiANC.freeze.10b.chr*.pass_only_phased.bcftools.vcf.gz | awk '{print $NF}'); do zgrep '^#' ${myfile} > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/$(basename ${myfile} ".vcf.gz").snpid.vcf; paste <(zgrep -v '^#' ${myfile} |awk '{print $1"\t"$2"\t"$1":"$2":"$4":"$5}') <(zgrep -v '^#' ${myfile} | cut -f4-) >> /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/$(basename ${myfile} ".vcf.gz").snpid.vcf; done &
## 
zgrep '^#' ${bcf_outfile}  > ${outDir}/$(basename ${bcf_outfile} ".vcf.gz").snpid.vcf

paste <(zgrep -v '^#' ${bcf_outfile} |awk '{print $1"\t"$2"\t"$1":"$2":"$4":"$5}') <(zgrep -v '^#' ${bcf_outfile} | cut -f4-) >> ${outDir}/$(basename ${bcf_outfile} ".vcf.gz").snpid.vcf

bgzip -f ${outDir}/$(basename ${bcf_outfile} ".vcf.gz").snpid.vcf

### Plink
## --set-all-var-ids '@:#:$1:$2'
## --new-id-max-allele-len 1000
i
## for chr in {10..22}; do /medpop/esp2/mesbah/tools/plink2_linux_x86_64_20200831 --vcf /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/topmed_N74k.MultiANC.freeze.10b.chr${chr}.pass_only_phased.bcftools.snpid.vcf.gz --maf 0.01 --hwe 1e-50 --geno 0.10 --mind 0.10 --threads $(nproc --all) --make-bed --out /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/plink/topmed_N74k.MultiANC.freeze.10b.chr${chr}.plink_maf01hwe1e50geno10mind10; done
###
${plink2} \
	--vcf ${outDir}/$(basename ${bcf_outfile} ".vcf.gz").snpid.vcf.gz \
	--maf 0.01 \
	--hwe 1e-50 \
	--geno 0.10 \
	--mind 0.10 \
	--threads $(nproc --all) \
	--make-bed \
	--out ${plink_outfile} 

###

###
rm ${bcf_outfile}
###

