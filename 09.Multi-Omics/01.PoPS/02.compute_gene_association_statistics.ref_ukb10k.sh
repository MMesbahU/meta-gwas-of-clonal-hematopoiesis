#!/bin/bash

## while read genename; do qsub -t 1-22 -R y -wd /broad/hptmp/mesbah/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=20G -l h_rt=30:00:00 -N ukbEUR_magma_ref_ukb10k.${genename} /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/PoPs/02.compute_gene_association_statistics.ref_ukb10k.sh /broad/hptmp/mesbah/ukb_chip/meta_gwas/pops/data/magma_0kb.genes.annot  /broad/hptmp/mesbah/ukb_chip/meta_gwas/pops/gwas/ukb/magma_in/ukb_eur.magma_summary.chr1_22.has${genename}.21Aug2021_ukbEUR.tsv "snp-wise=mean" /broad/hptmp/mesbah/ukb_chip/meta_gwas/pops/gwas/ukb/magma_out_refukb10k ${genename}; done < <(echo -e "CHIP\nDNMT3A\nTET2")

# bfile=${1}

gene_annot=${1}

pval=${2}

gene_model=${3}

out_dir=${4}

genename=${5}

chr=${SGE_TASK_ID}

out_prefix=${out_dir}/magma_out.ukb_eur.chr${chr}.has${genename}

bfile=/broad/hptmp/mesbah/ukb10k/ukb10k_eur.imp_chr${chr}_v3

######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
#########

## Run MAGMA
/medpop/esp2/mesbah/tools/MAGMA/magma \
	--bfile ${bfile} \
	--gene-annot ${gene_annot} \
	--pval ${pval} ncol=N \
	--gene-model ${gene_model} \
	--out ${out_prefix} \
	--batch ${chr} chr

######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

