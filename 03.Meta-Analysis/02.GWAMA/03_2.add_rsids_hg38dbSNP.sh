#!/bin/bash

## EUR: while read trait; do qsub -R y -wd /broad/hptmp/mesbah/gwas/ukb450k/metagwas/tmpdir -N add_rsids.EUR.${trait} -pe smp 1 -binding linear:1 -l h_rt=20:00:00 -l h_vmem=20G /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/03_2.add_rsids_hg38dbSNP.sh /medpop/esp/mesbah/GWAMA_VCF/hg38/eur_metaGWAS.${trait}.GWAMA.hg38.vcf.gz /broad/hptmp/mesbah/RefSeq/dbSNP/hg38/chr1_X.hg38.GCF_000001405.39.vcf.gz /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/eur_metaGWAS.${trait}.GWAMA.hg38_dbSNP.vcf.gz; done < <(echo -e "CHIP\nDNMT3A\nTET2")

# Run: while read trait; do qsub -R y -wd /broad/hptmp/mesbah/gwas/ukb450k/metagwas/tmpdir -N add_rsids.multiAncestry.${trait} -pe smp 1 -binding linear:1 -l h_rt=20:00:00 -l h_vmem=20G /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/03_2.add_rsids_hg38dbSNP.sh /medpop/esp/mesbah/GWAMA_VCF/hg38/metaGWAS.${trait}.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.hg38.vcf.gz /broad/hptmp/mesbah/RefSeq/dbSNP/hg38/chr1_X.hg38.GCF_000001405.39.vcf.gz /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/metaGWAS.${trait}.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.hg38_dbSNP.vcf.gz; done < <(echo -e "CHIP\nDNMT3A\nTET2")

###
source /broad/software/scripts/useuse

use Tabix

use Bcftools
 
######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
#########
annot_VCF=${1} # annotated GWAMA VCF

dbSNP_hg38=${2}

outVCFhg38=${3}
#################################
#### Add rsids from dbSNP 155
# bcftools annotate -a /broad/hptmp/mesbah/RefSeq/dbSNP/hg38/chr1_X.hg38.GCF_000001405.39.vcf.gz -c ID,GENEINFO,NSF,NSM,NSN,SYN,ASS,DSS,dbSNPBuildID,SSR /medpop/esp/mesbah/GWAMA_VCF/hg38/metaGWAS.CHIP.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.hg38.vcf.gz -Oz -o /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/metaGWAS.CHIP.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.hg38_dbSNP.v2.vcf.gz &
##
bcftools annotate \
	-a ${dbSNP_hg38} \
	-c ID,GENEINFO,NSF,NSM,NSN,SYN,ASS,DSS,dbSNPBuildID,SSR ${annot_VCF} \
	-Oz \
	-o ${outVCFhg38}
# 
tabix -p vcf ${outVCFhg38}

######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

