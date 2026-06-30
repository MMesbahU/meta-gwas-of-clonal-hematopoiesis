#!/bin/bash

### March 4, 2024
## Prep Samples to keep for filtering:
## while read ancestry; do  awk -F ' ' '{if(FILENAME==ARGV[1]){ids[$1]=$1}{if(FILENAME==ARGV[2] && ids[$2]){print $0}}}'  <(zcat /medpop/esp2/skoyama/passing/ukb_kgpprojection/v03/out/ukb.kgp_projected.tsv.gz| awk -v anc=${ancestry} 'NR>1 && $12==anc{print $1}') /medpop/esp2/pradeep/UKBiobank/v2data/ukb708_cal_chr1_v2_s488374.fam > /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/ukb708_cal_chr1_v2_s488374.knnGenAnc_${ancestry}.fam; done < <(zcat /medpop/esp2/skoyama/passing/ukb_kgpprojection/v03/out/ukb.kgp_projected.tsv.gz| awk 'NR>1{print $12}' | sort -V | uniq)
##
## Stratified by genetic ancestry calculated by SK. "knn" method
## while read ancestry; do qsub -R y -t 1-22 -wd /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/geno_qc/tmpdir  -N ${ancestry}.qc_ukb500k_geno -l h_vmem=20G -l h_rt=4:00:00 -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/0.prep4GWAS/scripts/qc_geno_ukb_v2.sh /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/geno_qc 1e-20 /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/ukb708_cal_chr1_v2_s488374.knnGenAnc_${ancestry}.fam ${ancestry}; done < <(echo -e "AFR\nAMR\nEAS\nEUR\nSAS")
## MultiAnc: while read ancestry; do qsub -R y -t 1-22 -wd /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/geno_qc/tmpdir  -N ${ancestry}.qc_ukb500k_geno -l h_vmem=20G -l h_rt=4:00:00 -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/0.prep4GWAS/scripts/qc_geno_ukb_v2.sh /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/geno_qc 1e-20 /medpop/esp2/pradeep/UKBiobank/v2data/ukb708_cal_chr1_v2_s488374.fam ${ancestry}; done < <(echo -e "MultiAnc") 

######################## 
## for chr in {1..22}; do qsub -wd /broad/hptmp/mesbah/gwas/ukb450k -N qc_ukb500k_geno_chr${chr} -l h_vmem=20G -l h_rt=1:00:00 -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/0.prep4GWAS/scripts/qc_geno_ukb_v2.sh ${chr} /broad/hptmp/mesbah/gwas/ukb450k; done
###
# source /broad/software/scripts/useuse
outDir=${1}
p_hwe=${2}
samples_2_keep=${3}
cohort=${4}
## 
chr=${SGE_TASK_ID}
out_file_name=${outDir}/qc_pass.${cohort}.chr${chr}
##
/medpop/esp2/mesbah/tools/plink2_linux_x86_64_20200831 \
	--bed /broad/ukbb/genotype/ukb_cal_chr${chr}_v2.bed \
	--bim /broad/ukbb/genotype/ukb_snp_chr${chr}_v2.bim \
	--fam /medpop/esp2/pradeep/UKBiobank/v2data/ukb708_cal_chr1_v2_s488374.fam \
	--maf 0.01 \
	--geno 0.1 \
	--hwe ${p_hwe} \
	--mind 0.1 \
	--keep ${samples_2_keep} \
	--write-snplist \
	--write-samples \
	--no-id-header \
	--out ${out_file_name}


############

### samples with >10 missing genotype in any autosome
## while read ANC; do cat /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/geno_qc/qc_pass.${ANC}.chr*.mindrem.id | sort -V | uniq > /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/geno_qc/qc_failed.${ANC}.chr1_22.mindrem.id; done < <(echo -e "AFR\nAMR\nSAS\nEAS\nEUR\nMultiAnc")

## samples to keep in GWAS null step: 
### '$1>0' is used to exclude negative eids in UKB
# while read ANC; do awk -F ' ' '{if(FILENAME==ARGV[1]){ids[$1]=$1}{if(FILENAME==ARGV[2] && !ids[$2]){print $0}}}' /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/geno_qc/qc_failed.${ANC}.chr1_22.mindrem.id /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/geno_qc/qc_pass.${ANC}.chr1.id | awk '$1>0' > /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/qc_pass_ukb_${ANC}_samples.v2.id; done < <(echo -e "AFR\nAMR\nSAS\nEAS\nEUR\nMultiAnc")

# or using grep: while read ANC; do grep -vf /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/geno_qc/qc_failed.${ANC}.chr1_22.mindrem.id /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/geno_qc/qc_pass.${ANC}.chr1.id | awk '$1>0' > /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/qc_pass_ukb_${ANC}_samples.id; done < <(echo -e "AFR\nAMR\nSAS\nEAS\nEUR\nMultiAnc")

## SNP to keep in null step:
# while read ANC; do for chr in {1..22}; do cat /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/geno_qc/qc_pass.${ANC}.chr${chr}.snplist  >> /medpop/esp/mesbah/GWAS_CHIP/inputUKB/2024_ukbb/geno_qc/qc_pass_snps.${ANC}.chr1_22.maf01_hwe1e20_geno10pct.snplist; done; done < <(echo -e "AFR\nAMR\nSAS\nEAS\nEUR\nMultiAnc")

############


## grep -vf <(cat /broad/hptmp/mesbah/gwas/ukb450k/qc_pass.chr*.mindrem.id | sort -V | uniq) /broad/hptmp/mesbah/gwas/ukb450k/qc_pass.chr1.id | awk '$1>0' >> /broad/hptmp/mesbah/gwas/ukb450k/qc_pass_ukb_all_samples.id

## for chr in {1..22}; do cat /broad/hptmp/mesbah/gwas/ukb450k/qc_pass.chr${chr}.snplist >> /broad/hptmp/mesbah/gwas/ukb450k/qc_pass_snps.ukb_v2.list; done

## >=10% genotype missing in any of the 22 autosomes
## 
## cat /broad/hptmp/mesbah/gwas/ukb450k/qc_pass.chr*.mindrem.id | sort -V | uniq > /broad/hptmp/mesbah/gwas/ukb450k/qc_failed.any_autosome_10.list

