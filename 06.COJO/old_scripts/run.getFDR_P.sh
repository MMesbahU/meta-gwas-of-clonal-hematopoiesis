#!/bin/bash

# qsub -wd /broad/hptmp/mesbah/ukb_chip/tmpdir -R y -N chip_meta4.getmy_fdr -l h_vmem=20G -l h_rt=02:00:00 -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.COJO/run.getFDR_P.sh /broad/hptmp/mesbah/ukb_chip/cojo/cojo_input/in_cojo_eaf001_min2Studies.8March2022.lifted_hg37.GWAMA.meta_CHIP.allsamples_topmed_ukbb_mgbb_BioVU.tsv
# qsub -wd /broad/hptmp/mesbah/ukb_chip/tmpdir -R y -N dnmt3a_meta4.getmy_fdr -l h_vmem=20G -l h_rt=02:00:00 -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.COJO/run.getFDR_P.sh /broad/hptmp/mesbah/ukb_chip/cojo/cojo_input/in_cojo_eaf001_min2Studies.8March2022.lifted_hg37.GWAMA.meta_DNMT3A.allsamples_topmed_ukbb_mgbb_BioVU.tsv

source /broad/software/scripts/useuse

use R-4.1


inputGWAS=${1}  # with column names SNP and p

Rscript /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.COJO/getFDR_P.R ${inputGWAS}


