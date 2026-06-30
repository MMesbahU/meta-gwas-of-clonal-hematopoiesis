#!/bin/bash

######################################### Combined Plink file ######################################
	## File plink list
# for chr in {2..22}; do echo "/broad/ukbb/genotype/ukb_cal_chr${chr}_v2.bed /broad/ukbb/genotype/ukb_snp_chr${chr}_v2.bim /medpop/esp2/pradeep/UKBiobank/v2data/ukb708_cal_chr1_v2_s488374.fam" >> /broad/hptmp/mesbah/gwas/ukb450k/ukb_list_beds.txt; done
for chr in {2..22}; do echo "/broad/ukbb/genotype/ukb_cal_chr${chr}_v2.bed /broad/ukbb/genotype/ukb_snp_chr${chr}_v2.bim /medpop/esp2/pradeep/UKBiobank/v2data/ukb708_cal_chr1_v2_s488374.fam" >> /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/ukb_list_beds.txt; done

	## Combine all autosomes
# qsub -R y -wd /broad/hptmp/mesbah/gwas/ukb450k -N combine_plink_beds -l h_vmem=30G -l h_rt=5:00:00 -pe smp 6 -binding linear:6  /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/0.prep4GWAS/scripts/combine_ukb_BEDs.sh 1 6 /broad/hptmp/mesbah/gwas/ukb450k/ukb_list_beds.txt /broad/hptmp/mesbah/gwas/ukb450k
qsub -R y -wd /broad/hptmp/mesbah/dataset/ukb450 -N combine_plink_beds -l h_vmem=30G -l h_rt=5:00:00 -pe smp 6 -binding linear:6  /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/0.prep4GWAS/scripts/combine_ukb_BEDs.sh 1 6 /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/ukb_list_beds.txt /broad/hptmp/mesbah/dataset/ukb450
########################################################################################################


######################################## Genotype QC for Regenie Step 1 #########################
	## Basic QCs per Autosome:
	for chr in {1..22}; do qsub -wd /broad/hptmp/mesbah/gwas/ukb450k -N qc_ukb500k_geno_chr${chr} -l h_vmem=20G -l h_rt=1:00:00 -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/0.prep4GWAS/scripts/qc_geno_ukb_v2.sh ${chr} /broad/hptmp/mesbah/gwas/ukb450k; done

	## qc failed Samples
	## >=10% genotype missing in any of the 22 autosomes
        cat /broad/hptmp/mesbah/gwas/ukb450k/qc_pass.chr*.mindrem.id | sort -V | uniq > /broad/hptmp/mesbah/gwas/ukb450k/qc_failed.any_autosome_10.list
	## qc passed Samples
	grep -vf /broad/hptmp/mesbah/gwas/ukb450k/qc_failed.any_autosome_10.list /broad/hptmp/mesbah/gwas/ukb450k/qc_pass.chr1.id | awk '$1>0' >> /broad/hptmp/mesbah/gwas/ukb450k/qc_pass_ukb_all_samples.id
	## qc passed SNPs
	for chr in {1..22}; do cat /broad/hptmp/mesbah/gwas/ukb450k/qc_pass.chr${chr}.snplist >> /broad/hptmp/mesbah/gwas/ukb450k/qc_pass_snps.ukb_v2.list; done

	## >=10% genotype missing in any of the 22 autosomes
	cat /broad/hptmp/mesbah/gwas/ukb450k/qc_pass.chr*.mindrem.id | sort -V | uniq > /broad/hptmp/mesbah/gwas/ukb450k/qc_failed.any_autosome_10.list

###################################################################################################	



