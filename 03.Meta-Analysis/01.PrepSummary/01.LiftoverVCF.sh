#!/bin/bash

# qsub -R y -wd /broad/hptmp/mesbah/tmpdir -N Liftover.chip -pe smp 1 -binding linear:1 -l h_rt=20:00:00 -l h_vmem=20G /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/01.PrepSummary/01.LiftoverVCF.sh /medpop/esp2/mesbah/tools/LiftOVer/picard.jar /broad/hptmp/mesbah/ukb_chip/eurUKB/ukb_Chr1_22.hasCHIP.regenie.vcf.gz /broad/hptmp/mesbah/ukb_chip/eurUKB/lifted.ukb_Chr1_22.hasCHIP.regenie.vcf.gz /medpop/esp2/mesbah/tools/LiftOVer/b37ToHg38.over.chain /broad/hptmp/mesbah/ukb_chip/eurUKB/rejected_variants.ukb_Chr1_22.hasCHIP.regenie.vcf.gz /medpop/esp/skoyama/publicdata/reference/Homo_sapiens_assembly38.fasta 

source /broad/software/scripts/useuse

use Picard-Tools

### input files
picard=${1}
inVCF=${2}
outVCF=${3}
chain=${4}
rejectVCF=${5}
refTarget=${6}

## Prepare input VCF


## Run Picard to liftover
java -jar ${picard} LiftoverVcf \
	-I ${inVCF} \
	-O ${outVCF} \
	-C ${chain} \
	--REJECT ${rejectVCF} \
	-R ${refTarget}

## Prepare Output summary File for meta-analysis


