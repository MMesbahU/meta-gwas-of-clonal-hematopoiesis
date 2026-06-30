#!/bin/bash

#
## 2024 March 4:
# qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/tmpdir -N combine_plink_beds -l h_vmem=20G -l h_rt=5:00:00 -pe smp 2 -binding linear:2 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/0.prep4GWAS/topmed/02.merge_plink_files.sh 1 2 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/plink/bed_list.txt /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/plink
##


source /broad/software/scripts/useuse

use .plink-1.90b

# input
chr=${1} # 1

numThreads=${2}

list_beds=${3} #

outDir=${4}

plinkBEDfile=${outDir}/topmed_N74k.MultiANC.freeze.10b.chr${chr}.plink_maf01hwe1e50geno10mind10.bed

plinkBIMfile=${outDir}/topmed_N74k.MultiANC.freeze.10b.chr${chr}.plink_maf01hwe1e50geno10mind10.bim

plinkFAMfile=${outDir}/topmed_N74k.MultiANC.freeze.10b.chr${chr}.plink_maf01hwe1e50geno10mind10.fam

outPlinkfile=${outDir}/topmed_N74k.MultiANC.freeze.10b.chr1_22.plink_maf01hwe1e50geno10mind10
###### 
rm ${list_beds}
for chrom in {2..22}; do echo "/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/plink/topmed_N74k.MultiANC.freeze.10b.chr${chrom}.plink_maf01hwe1e50geno10mind10.bed /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/plink/topmed_N74k.MultiANC.freeze.10b.chr${chrom}.plink_maf01hwe1e50geno10mind10.bim /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/geno/plink/topmed_N74k.MultiANC.freeze.10b.chr${chrom}.plink_maf01hwe1e50geno10mind10.fam" >> ${list_beds}; done

############################# Clock Time: Start
echo -e "Job started at: $(date)"

Job_START=$(date +%s)
############################

#plink2=/medpop/esp2/mesbah/tools/plink2
# ${plink2} 
plink \
	--bed ${plinkBEDfile} \
	--bim ${plinkBIMfile} \
	--fam ${plinkFAMfile} \
	--merge-list ${list_beds} \
	--threads $(nproc --all) \
	--make-bed \
	--out ${outPlinkfile}

############################# Clock Time: END
echo "Job ended at: $(date)" 

Job_END=$(date +%s)


echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'

##############################


