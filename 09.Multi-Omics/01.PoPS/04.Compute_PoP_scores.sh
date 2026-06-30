#!/bin/bash

source /broad/software/scripts/useuse

reuse Anaconda3

### 24 May, 2022
### POPS with all meta-GWAS
## while read genename; do qsub -R y -t 1-22 -wd /broad/hptmp/mesbah/ukb_chip/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=60G -l h_rt=30:00:00 -N step4.POPs.${genename}.refUKB200k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/PoPs/04.Compute_PoP_scores.sh /medpop/esp2/projects/software/PoPS/data/PoPS.features.txt.gz /medpop/esp2/projects/software/PoPS/data/gene_loc.txt /broad/hptmp/mesbah/ukb_chip/PoPs/output_pops/combined_output/chr1_22.lifted_hg37.GWAMA.meta_${genename}.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies /broad/hptmp/mesbah/ukb_chip/PoPs/output_pops/combined_output/chr1_22.pops_step3.has${genename}.features /medpop/esp2/projects/software/PoPS/data/control.features /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/PoPs/pops/pops.predict_scores.py /broad/hptmp/mesbah/ukb_chip/PoPs/output_pops/combined_output/final_pops_step4.has${genename}; done < <(echo -e "CHIP\nDNMT3A\nTET2")

########

## 30 Dec, 2021
## while read genename; do qsub -R y -t 1-22 -wd /broad/hptmp/mesbah/ukb_chip/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=60G -l h_rt=30:00:00 -N step4.ukb200k.${genename}.ukb200kREF /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/PoPs/04.Compute_PoP_scores.sh /broad/hptmp/mesbah/ukb_chip/PoPs/pops_features/PoPS.features.txt.gz /broad/hptmp/mesbah/ukb_chip/PoPs/pops_features/gene_loc.txt /broad/hptmp/mesbah/ukb_chip/PoPs/output_pops/chr1_22.magma_out.ukb_all_200k.has${genename} /broad/hptmp/mesbah/ukb_chip/PoPs/output_pops/chr1_22.pops_step3.has${genename}.features /broad/hptmp/mesbah/ukb_chip/PoPs/pops_features/control.features /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/PoPs/pops/pops.predict_scores.py /broad/hptmp/mesbah/ukb_chip/PoPs/output_pops/final_pops_step4.has${genename}; done < <(echo -e "CHIP\nDNMT3A\nTET2")

## while read genename; do qsub -R y -t 1-22 -wd /broad/hptmp/mesbah/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=60G -l h_rt=30:00:00 -N step3.ukbEUR_magma.${genename}.ref1kg /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/PoPs/04.Compute_PoP_scores.sh /broad/hptmp/mesbah/ukb_chip/meta_gwas/pops/data/PoPS.features.txt.gz /broad/hptmp/mesbah/ukb_chip/meta_gwas/pops/data/gene_loc.txt /broad/hptmp/mesbah/ukb_chip/meta_gwas/pops/gwas/ukb/magma_out/magma_out.ukb_eur.chr1_22.has${genename} /broad/hptmp/mesbah/ukb_chip/meta_gwas/pops/gwas/ukb/magma_out/pops_step2.has${genename}.features /broad/hptmp/mesbah/ukb_chip/meta_gwas/pops/data/control.features /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/PoPs/pops/pops.predict_scores.py /broad/hptmp/mesbah/ukb_chip/meta_gwas/pops/gwas/ukb/magma_out/pops_step3.has${genename}; done < <(echo -e "CHIP\nDNMT3A\nTET2")

##
pops_features=${1}

gene_loc=${2}

step1_gene_results=${3}

step2_selected_features=${4}

control_features=${5}

pops_scripts=${6} # /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/PoPs/pops/pops.predict_scores.py

pops_output=${7}

chr=${SGE_TASK_ID}


# out_prefix=${out_dir}/magma_out.ukb_eur.chr${chr}.has${genename}

######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
#########

## Run MAGMA
python ${pops_scripts} \
	--gene_loc ${gene_loc} \
	--gene_results ${step1_gene_results} \
	--features ${pops_features} \
	--selected_features ${step2_selected_features} \
	--control_features ${control_features} \
	--chromosome ${chr} \
	--out ${pops_output}

######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

