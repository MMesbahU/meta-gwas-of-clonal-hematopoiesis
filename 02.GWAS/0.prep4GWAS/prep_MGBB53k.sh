#!/bin/bash

	# Directly Genotyped QC data: chr1-22
	## ANC Stratified
# while read ANC; do /medpop/esp2/mesbah/tools/plink2_linux_x86_64_20200831 --bfile  /medpop/esp2/projects/MGB_Biobank/genotype/53K_GSA/release/GSA_53K.hg38 --exclude <(awk '$1>22 {print $2}' /medpop/esp2/projects/MGB_Biobank/genotype/53K_GSA/release/GSA_53K.hg38.bim) --keep <( grep -f <( zcat /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.${ANC}.Feb2024.tsv.gz | awk '{print $3}') /medpop/esp2/projects/MGB_Biobank/genotype/53K_GSA/release/GSA_53K.hg38.fam) --make-bed --maf 0.01 --hwe 1e-15 --geno 0.1 --out /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/qced.mgbb53k.autosomes.${ANC}.hew1e15maf01; done < <(echo -e "AFR\nSAS\nEAS") 1>>/medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/prep.MGB_geno.log 2>>/medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/prep.MGB_geno.err &

# while read ANC; do /medpop/esp2/mesbah/tools/plink2_linux_x86_64_20200831 --bfile  /medpop/esp2/projects/MGB_Biobank/genotype/53K_GSA/release/GSA_53K.hg38 --exclude <(awk '$1>22 {print $2}' /medpop/esp2/projects/MGB_Biobank/genotype/53K_GSA/release/GSA_53K.hg38.bim) --keep <( grep -f <( zcat /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.${ANC}.Feb2024.tsv.gz | awk '{print $3}') /medpop/esp2/projects/MGB_Biobank/genotype/53K_GSA/release/GSA_53K.hg38.fam) --write-snplist --write-samples --no-id-header --maf 0.01 --hwe 1e-15 --geno 0.1 --out /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/qced.mgbb53k.autosomes.${ANC}.hew1e15maf01; done < <(echo -e "AFR\nSAS\nEAS") 1>>/medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/prep.MGB_geno.log 2>>/medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/prep.MGB_geno.err &

	# Feb 2024
	/medpop/esp2/mesbah/tools/plink2_linux_x86_64_20200831 --bfile  /medpop/esp2/projects/MGB_Biobank/genotype/53K_GSA/release/GSA_53K.hg38 --exclude <(awk '$1>22 {print $2}' /medpop/esp2/projects/MGB_Biobank/genotype/53K_GSA/release/GSA_53K.hg38.bim) --make-bed --maf 0.01 --mac 100 --hwe 1e-20 --geno 0.1 --out /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/qced.mgbb53k.autosomes.hew1e20maf01mac100 1>>/medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/prep.MGB_geno.log 2>>/medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/prep.MGB_geno.err &
	## 53,306 samples, 417,424 variants passed qc

	# old: /medpop/esp2/mesbah/tools/plink2_linux_x86_64_20200831 --bfile  /medpop/esp2/projects/MGB_Biobank/genotype/53K_GSA/release/GSA_53K.hg38 --exclude <(awk '$1>22 {print $2}' /medpop/esp2/projects/MGB_Biobank/genotype/53K_GSA/release/GSA_53K.hg38.bim) --make-bed --maf 0.01 --mac 100 --hwe 1e-15 --geno 0.1 --out /broad/hptmp/mesbah/gwas/mgbb53k/qced.mgbb53k.autosomes
	# 53,306 Samples and 404,770 SNPs passed qc

	## PCA 
	qsub -b y -R y -wd /broad/hptmp/mesbah/gwas/mgbb53k -pe smp 2 -binding linear:2 -l h_vmem=20G -l h_rt=10:00:00 -N prep_mgb53k_chr1_22 /medpop/esp2/mesbah/tools/plink2_linux_x86_64_20200831 --bfile /broad/hptmp/mesbah/gwas/mgbb53k/qced.mgbb53k.autosomes --pca 20 --out /broad/hptmp/mesbah/gwas/mgbb53k/PC1_20.qced.mgbb53k.autosomes --threads 2

	## KING: did not finish  
	# 3rd degree
	## https://www.kingrelatedness.com/manual.shtml: "Please do not prune or filter any "good" SNPs that pass QC prior to any KING inference, unless the number of variants is too many to fit the computer memory, e.g., > 100,000,000 as in a WGS study, in which case rare variants can be filtered out. LD pruning is not recommended in KING."
	# qsub -b y -R y -wd /broad/hptmp/mesbah/gwas/mgbb53k -pe smp 8 -binding linear:8 -l h_vmem=5G -l h_rt=40:00:00 -N king_mgbb53k /medpop/esp2/mesbah/tools/king -b /broad/hptmp/mesbah/gwas/mgbb53k/qced.mgbb53k.autosomes.bed --kinship --degree 3 --cpus 8 --prefix /broad/hptmp/mesbah/gwas/mgbb53k/KIN3rd.qced.mgbb53k.autosomes
	qsub -R y -wd /broad/hptmp/mesbah/gwas/mgbb53k -pe smp 8 -binding linear:8 -l h_vmem=5G -l h_rt=40:00:00 -N king_mgbb53k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/0.prep4GWAS/scripts/king.sh /broad/hptmp/mesbah/gwas/mgbb53k/qced.mgbb53k.autosomes.bed /broad/hptmp/mesbah/gwas/mgbb53k/KIN3rd.qced.mgbb53k.autosomes 8
	
	## WHITE:
	awk '(NR==1){print $0}(NR>1 && $14=="WHITE"){print $0}' /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.21Jul2022.tsv > /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k_WHITE.imp_new.noRel_sk.21Jul2022.tsv

