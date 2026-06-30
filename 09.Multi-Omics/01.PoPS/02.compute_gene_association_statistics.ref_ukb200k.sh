#!/bin/bash

##### Aug 26, 2022
############################# 650k meta-GWAS
# /medpop/esp2/mesbah/tools/MAGMA/magma --bfile /broad/hptmp/mesbah/ukb_chip/ukb200k_reference/plink_ukb200k/ukb_imp_chr22_v3.ukb200k_CHIP_samples --gene-annot /medpop/esp2/projects/software/PoPS/data/magma_0kb.genes.annot --pval /broad/hptmp/mesbah/gwas/cojo/input/ref_ukb200k.summary_4_cojo.lifted_hg37.metaGWAS.CHIP.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.hg37_dbSNP.tsv ncol=N --gene-model snp-wise=mean --batch 22 chr --out test.chr22.magma.multi-CHIP

## Multi-ancestry:
# mkdir -p /broad/hptmp/mesbah/gwas/pops/{tmpdir,outdir}
# while read genename; do qsub -t 1-22 -R y -wd /broad/hptmp/mesbah/gwas/pops/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=20G -l h_rt=20:00:00 -N step1.${genename}.magma.multiAncestry.ukb200k_Ref /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/10.Multi-Omics/01.PoPS/02.compute_gene_association_statistics.ref_ukb200k.sh /medpop/esp2/projects/software/PoPS/data/magma_0kb.genes.annot /broad/hptmp/mesbah/gwas/cojo/input/ref_ukb200k.summary_4_cojo.lifted_hg37.metaGWAS.${genename}.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.hg37_dbSNP.tsv "snp-wise=mean" /broad/hptmp/mesbah/gwas/pops/outdir ${genename} /broad/hptmp/mesbah/ukb_chip/ukb200k_reference/plink_ukb200k; done < <(echo -e "CHIP\nDNMT3A\nTET2")

## EUR only GWAS:
# while read genename; do qsub -t 1-22 -R y -wd /broad/hptmp/mesbah/gwas/pops/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=20G -l h_rt=20:00:00 -N step1.${genename}.magma.eur.ukb200k_Ref /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/10.Multi-Omics/01.PoPS/02.compute_gene_association_statistics.ref_ukb200k.sh /medpop/esp2/projects/software/PoPS/data/magma_0kb.genes.annot /broad/hptmp/mesbah/gwas/cojo/input/ref_ukb200k.summary_4_cojo.lifted_hg37.eur_metaGWAS.${genename}.GWAMA.hg37_dbSNP.tsv "snp-wise=mean" /broad/hptmp/mesbah/gwas/pops/outdir ${genename} /broad/hptmp/mesbah/ukb_chip/ukb200k_reference/plink_ukb200k; done < <(echo -e "CHIP\nDNMT3A\nTET2")

##############################################


## Copy PoPs features
## cp /medpop/esp2/akhil/topmed/lipid_mashr/post_betas/{PoPS.features.txt.gz,control.features,gene_loc.txt,magma_0kb.genes.annot} /broad/hptmp/mesbah/ukb_chip.v2/PoPs/pops_features/

######################## Meta-GWAS all TOPMed2019 + UKB200k + MGBB 13k + BioVU 54k
### 18 May 2022
### while read genename; do qsub -t 1-22 -R y -wd /broad/hptmp/mesbah/ukb_chip/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=20G -l h_rt=20:00:00 -N ${genename}.meta4.n2maf001.ukb200k_Ref /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/PoPs/02.compute_gene_association_statistics.ref_ukb200k.sh /medpop/esp2/projects/software/PoPS/data/magma_0kb.genes.annot /broad/hptmp/mesbah/ukb_chip/PoPs/PoPs_input/chr1_22.lifted_hg37.GWAMA.meta_${genename}.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.tsv "snp-wise=mean" /broad/hptmp/mesbah/ukb_chip/PoPs/output_pops ${genename} /broad/hptmp/mesbah/ukb_chip/ukb200k_reference/plink_ukb200k; done < <(echo -e "CHIP\nDNMT3A\nTET2")
####################################################
###########

## Jan 20, 2022
### from COJO input
## while read gene; do echo -e "SNP\tP\tN" > PoPs/PoPs_input/${gene}.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.tsv; done < <(echo -e "CHIP\nDNMT3A\nTET2")
## while read gene; do awk 'NR>1{print $1"\t"$7"\t"$8}' cojo/cojo_input/lifted_hg37.GWAMA.meta_${gene}.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.cojo.tsv >> PoPs/PoPs_input/${gene}.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.tsv; done < <(echo -e "CHIP\nDNMT3A\nTET2") &

## while read genename; do qsub -t 1-22 -R y -wd /broad/hptmp/mesbah/ukb_chip.v2/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=20G -l h_rt=20:00:00 -N ALL_magma_ref_ukb200k.${genename} /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/PoPs/02.compute_gene_association_statistics.ref_ukb200k.sh /broad/hptmp/mesbah/ukb_chip.v2/PoPs/pops_features/magma_0kb.genes.annot /broad/hptmp/mesbah/ukb_chip.v2/PoPs/PoPs_input/${genename}.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.tsv "snp-wise=mean" /broad/hptmp/mesbah/ukb_chip.v2/PoPs/output_pops ${genename}; done < <(echo -e "CHIP\nDNMT3A\nTET2")
############################

################### UKB200k GWAS only
## while read genename; do echo -e "SNP\tP\tN" > /broad/hptmp/mesbah/ukb_chip/PoPs/PoPs_input/ukb_all200k.${genename}.summary.tsv; done < <(echo -e "CHIP\nDNMT3A\nTET2");
# while read genename; do awk 'NR>1{print $1"\t"$7"\t"$8}' /broad/hptmp/mesbah/ukb_chip/cojo/cojo_input/summary_for_cojo.maf0.1_info0.3.chr1_22.has${genename}.11Aug2021_ukb200k.tsv >> /broad/hptmp/mesbah/ukb_chip/PoPs/PoPs_input/ukb_all200k.${genename}.summary.tsv; done < <(echo -e "CHIP\nDNMT3A\nTET2") &

## while read genename; do qsub -t 1-22 -R y -wd /broad/hptmp/mesbah/ukb_chip/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=20G -l h_rt=20:00:00 -N ukbALL_magma_ref_ukb200k.${genename} /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/PoPs/02.compute_gene_association_statistics.ref_ukb200k.sh /broad/hptmp/mesbah/ukb_chip/PoPs/pops_features/magma_0kb.genes.annot /broad/hptmp/mesbah/ukb_chip/PoPs/PoPs_input/ukb_all200k.${genename}.summary.tsv "snp-wise=mean" /broad/hptmp/mesbah/ukb_chip/PoPs/output_pops ${genename}; done < <(echo -e "CHIP\nDNMT3A\nTET2")
#################

# bfile=${1}

gene_annot=${1}

pval=${2} ## file with SNP, P and N information

gene_model=${3}

out_dir=${4}

genename=${5}

chr=${SGE_TASK_ID}

out_prefix=${out_dir}/pops_step1_magma.chr${chr}.$(basename ${pval} ".tsv")

refDir=${6}

bfile=${refDir}/ukb_imp_chr${chr}_v3.ukb200k_CHIP_samples #/broad/hptmp/mesbah/ukb10k/ukb10k_eur.imp_chr${chr}_v3

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

