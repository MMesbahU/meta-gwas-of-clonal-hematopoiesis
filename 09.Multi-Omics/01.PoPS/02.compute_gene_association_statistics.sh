#!/bin/bash

## while read genename; do qsub -R y -wd /broad/hptmp/mesbah/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=20G -l h_rt=30:00:00 -N ukbEUR_magma_ref1k.${genename} /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/PoPs/02.compute_gene_association_statistics.sh /broad/hptmp/mesbah/ukb_chip/meta_gwas/pops/data/1000G.EUR /broad/hptmp/mesbah/ukb_chip/meta_gwas/pops/data/magma_0kb.genes.annot  /broad/hptmp/mesbah/ukb_chip/meta_gwas/pops/gwas/ukb/magma_in/ukb_eur.magma_summary.chr1_22.has${genename}.21Aug2021_ukbEUR.tsv "snp-wise=mean" /broad/hptmp/mesbah/ukb_chip/meta_gwas/pops/gwas/ukb/magma_out/magma_out.ukb_eur.chr1_22.has${genename}; done < <(echo -e "CHIP\nDNMT3A\nTET2")

bfile=${1}

gene_annot=${2}

pval=${3}

gene_model=${4}

out_prefix=${5}

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
	--out ${out_prefix}

# --batch 1 chr
######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

