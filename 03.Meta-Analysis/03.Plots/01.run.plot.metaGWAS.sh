#!/bin/bash

# N study>=3
# while read trait; do qsub -R y -wd /broad/hptmp/mesbah/dataset/ch_gwas/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=20G -l h_rt=10:00:00 -N plot_${trait}_700k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/01.run.plot.metaGWAS.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/01.plot.metaGWAS.R /broad/hptmp/mesbah/dataset/ch_gwas/for_plot.metaGWAS.has${trait}.topmed74k_topImpukb200k250k_aou98k_mgbb53k_biovu54k.nStd3pls.tsv.gz /broad/hptmp/mesbah/dataset/ch_gwas ${trait} maf001_nStd3pls; done < <(echo -e "CH\nDNMT3A\nTET2")
# N study>=2
# while read trait; do qsub -R y -wd /broad/hptmp/mesbah/dataset/ch_gwas/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=20G -l h_rt=10:00:00 -N plot_${trait}_700k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/01.run.plot.metaGWAS.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/01.plot.metaGWAS.R /broad/hptmp/mesbah/dataset/ch_gwas/for_plot.metaGWAS.has${trait}.topmed74k_topImpukb200k250k_aou98k_mgbb53k_biovu54k.tsv.gz /broad/hptmp/mesbah/dataset/ch_gwas ${trait}; done < <(echo -e "DNMT3A\nTET2")

#############
## Run: qsub -R y -wd /broad/hptmp/mesbah/ukb_chip/meta_gwas/eurGWAS_meta4/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=20G -l h_rt=10:00:00 -N plot_eur_meta4.chip.23Aug2021 /broad/hptmp/mesbah/ukb_chip/meta_gwas/eurGWAS_meta4/scripts/run.plot.metaGWAS.sh /broad/hptmp/mesbah/ukb_chip/meta_gwas/eurGWAS_meta4/scripts/plot.metaGWAS.R /broad/hptmp/mesbah/ukb_chip/meta_gwas/eurGWAS_meta4/locusZoom/eurMeta4_hasCHIP.ivw_Metal_random.TOPMed29373_UKB167713_MGB9950_BioVUsaige52640_1583.locusZoom_maf001.tsv.gz /broad/hptmp/mesbah/ukb_chip/meta_gwas/eurGWAS_meta4/locusZoom CHIP

## Run : qsub -R y -wd /broad/hptmp/mesbah/ukb_chip/meta_gwas/meta4.topukbmgbbiovuSaige/locusZoom -pe smp 1 -binding linear:1 -l h_vmem=20G -l h_rt=20:00:00 -N plot_metagwas.21Aug2021 /broad/hptmp/mesbah/ukb_chip/meta_gwas/meta4.topukbmgbbiovuSaige/scripts/run.plot.metaGWAS.sh /broad/hptmp/mesbah/ukb_chip/meta_gwas/meta4.topukbmgbbiovuSaige/scripts/plot.metaGWAS.R

source /broad/software/scripts/useuse

use R-4.0

plot_script=${1}

gwas_summary=${2}

Outdir=${3}

myTrait=${4}

pic_prefix=${5}

# Rscript ${plot_script} #/broad/hptmp/mesbah/ukb_chip/meta_gwas/meta4.topukbmgbbiovuSaige/scripts/plot.metaGWAS.R


Rscript ${plot_script} ${gwas_summary} ${Outdir} ${myTrait} ${pic_prefix}


