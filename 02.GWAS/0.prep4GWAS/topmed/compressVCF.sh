#!/bin/bash

##
## for chr in {1,10,11}; do qsub -r y -N chr${chr}.compress -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/tmpdir -l h_rt=20:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/0.prep4GWAS/topmed/compressVCF.sh ${chr}; done

source /broad/software/scripts/useuse

use Tabix

#
chr=${1}
# for chr in {1,10,11}; do qsub -b y -r y -N chr${chr}.compress -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/tmpdir gzip -f /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/topmed_N74k.MultiANC.freeze.10b.chr${chr}.pass_only_phased.bcftools.snpid.vcf; done 

## 
bgzip -f /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/topmed_N74k.MultiANC.freeze.10b.chr${chr}.pass_only_phased.bcftools.snpid.vcf

##
