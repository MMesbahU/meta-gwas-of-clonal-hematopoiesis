#!/bin/bash

#### Note: ## SORT step in LiftoverVcf requires lots of Memory!!!! 
# solution: DISABLE_SORT=true
#####

### Prepare Reference File
# wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.{fasta.gz,fasta.fai}
## Get reference from MOCHA
# wget https://software.broadinstitute.org/software/mocha/mocha.GRCh37.zip
# unzip mocha.GRCh37.zip human_g1k_v37.{fasta,fasta.fai}
# prepare Dict file
# use Picard-Tools; java -jar /medpop/esp2/mesbah/tools/LiftOVer/picard.jar CreateSequenceDictionary R=/broad/hptmp/mesbah/RefSeq/human_g1k_v37.fasta O=/broad/hptmp/mesbah/RefSeq/human_g1k_v37.dict
###################################
## Multi-Ancestry: 
# while read gene_name; do qsub -R y -wd /broad/hptmp/mesbah/gwas/lifted_hg37 -N lift_hg37_all.${gene_name} -pe smp 1 -binding linear:1 -l h_rt=20:00:00 -l h_vmem=20G /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/04.GWAMA_hg38_to_hg37.sh /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/metaGWAS.${gene_name}.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.hg38_dbSNP.vcf.gz /broad/hptmp/mesbah/gwas/lifted_hg37/lifted_hg37.metaGWAS.${gene_name}.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.hg37_dbSNP.vcf.gz /medpop/esp2/mesbah/tools/LiftOVer/grch38_to_grch37.over.chain /medpop/esp2/mesbah/tools/LiftOVer/picard.jar /broad/hptmp/mesbah/RefSeq/human_g1k_v37.fasta /broad/hptmp/mesbah/gwas/lifted_hg37/rejected_hg37.metaGWAS.${gene_name}.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.hg37_dbSNP.vcf.gz; done < <(echo -e "CHIP\nDNMT3A\nTET2")

## EUR only: 
# while read gene_name; do qsub -R y -wd /medpop/esp/mesbah/GWAMA_VCF/hg37 -N lift_hg37.eur.${gene_name} -pe smp 1 -binding linear:1 -l h_rt=20:00:00 -l h_vmem=20G /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/04.GWAMA_hg38_to_hg37.sh /medpop/esp/mesbah/GWAMA_VCF/hg38/eur_metaGWAS.${gene_name}.GWAMA.hg38_dbSNP.vcf.gz /medpop/esp/mesbah/GWAMA_VCF/hg37/lifted_hg37.eur_metaGWAS.${gene_name}.GWAMA.hg37_dbSNP.vcf.gz /medpop/esp2/mesbah/tools/LiftOVer/grch38_to_grch37.over.chain /medpop/esp2/mesbah/tools/LiftOVer/picard.jar /broad/hptmp/mesbah/RefSeq/human_g1k_v37.fasta /medpop/esp/mesbah/GWAMA_VCF/hg37/rejected_hg37.eur_metaGWAS.${gene_name}.GWAMA.hg37_dbSNP.vcf.gz; done < <(echo -e "CHIP\nDNMT3A\nTET2")

###############################
source /broad/software/scripts/useuse

use Picard-Tools

use Tabix

use Bcftools

use VCFtools
######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
#########
inVCF=${1}
outVCF=${2}
hg38tohg37=${3} # /medpop/esp2/mesbah/tools/LiftOVer/grch38_to_grch37.over.chain
Picard=${4} # /medpop/esp2/mesbah/tools/LiftOVer/picard.jar
refTarget=${5} #/broad/hptmp/mesbah/RefSeq/human_g1k_v37.fasta
rejectVCF=${6}
#
# java -Xmx20G -jar /medpop/esp2/mesbah/tools/LiftOVer/picard.jar LiftoverVcf I=/broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/metaGWAS.CHIP.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.hg38.vcf.gz O=/broad/hptmp/mesbah/gwas/lifted_hg37/v2.lifted_hg37.GWAMA_CHIP.vcf.gz REJECT=/broad/hptmp/mesbah/gwas/lifted_hg37/v2.rejected_hg37.GWAMA_CHIP.vcf.gz CHAIN=/medpop/esp2/mesbah/tools/LiftOVer/grch38_to_grch37.over.chain R=/broad/hptmp/mesbah/RefSeq/human_g1k_v37.fasta VERBOSITY=DEBUG DISABLE_SORT=true

java -Xmx20G -jar ${Picard} LiftoverVcf \
	I=${inVCF} \
	O=${outVCF} \
	CHAIN=${hg38tohg37} \
	REJECT=${rejectVCF} \
	R=${refTarget} \
	VERBOSITY=DEBUG \
	DISABLE_SORT=true

# Need sorting to index!!
# tabix -f -p vcf ${outVCF}

######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

