#!/bin/bash

## Get variants
for chr in {1..22}; do /medpop/esp2/mesbah/tools/bgenix/build/apps/bgenix -g /broad/ukbb/imputed_v3/ukb_imp_chr${chr}_v3.bgen -list > /medpop/esp2/mesbah/datasets/ukbb/ukbb_v3/var/hg37_var.ukb_imp_chr${chr}_v3.txt; done &

## hapmap v3 imputed snps for afr,gbr reference
# cat /broad/hptmp/mesbah/dataset/ch_gwas/ldak/ref/*.hapmap/*.hapmap.cors.bim | sort -V -k2 | uniq > /medpop/esp2/mesbah/datasets/ukbb/ukbb_v3/var/chr1_22.afr_or_gbr.ldak.var.hapmap.cors.bim

## Get overlapping snps
# for chr in {5..13}; do awk -v OFS='\t' 'NR==FNR { vars[$1]; next } FNR==1 || $1 in vars { print $2, $3, $4, $5, $6, $7, $8} ' <(awk -v mychr=${chr} '$1==mychr {print $2":"$5":"$6}' /medpop/esp2/mesbah/datasets/ukbb/ukbb_v3/var/chr1_22.afr_or_gbr.ldak.var.hapmap.cors.bim) <(awk -v mychr=${chr} 'NR>1{print mychr":"$4":"$6":"$7"\t"$0}' /medpop/esp2/mesbah/datasets/ukbb/ukbb_v3/var/hg37_var.ukb_imp_chr${chr}_v3.txt) > /medpop/esp2/mesbah/datasets/ukbb/ukbb_v3/var/hapmapV3_hg37_varid.ukb_imp_chr${chr}_v3.txt; gzip -f /medpop/esp2/mesbah/datasets/ukbb/ukbb_v3/var/hg37_var.ukb_imp_chr${chr}_v3.txt;done &

##
# head -n2 /broad/ukbb/imputed_v3/ukb_mfi_chr1_v3.txt
# 1:10177_A_AC	rs367896724	10177	A	AC	0.40079	AC	0.467935
# 1:10235_T_TA	rs540431307	10235	T	TA	0.000367353	TA	0.214688

for chr in {1..22}; do awk -v OFS='\t' 'NR==FNR { vars[$1]; next } $1 in vars { print $2, $3, $4, $5, $6, $7, $8} ' <(awk -v mychr=${chr} '$1==mychr {print $2":"$5":"$6}' /medpop/esp2/mesbah/datasets/ukbb/ukbb_v3/var/chr1_22.afr_or_gbr.ldak.var.hapmap.cors.bim) <(awk -v mychr=${chr} '{print mychr":"$3":"$4":"$5"\t"$0}' /broad/ukbb/imputed_v3/ukb_mfi_chr${chr}_v3.txt) > /medpop/esp2/mesbah/datasets/ukbb/ukbb_v3/var/hapmap3.ukb_mfi_chr${chr}_v3.txt; done

## GT




/medpop/esp2/mesbah/tools/bgenix/build/apps/bgenix -g /broad/ukbb/imputed_v3/ukb_imp_chr1_v3.bgen -vcf -incl-rsids /medpop/esp2/mesbah/datasets/ukbb/ukbb_v3/var/hapmapV3_hg37_varid.ukb_imp_chr1_v3.txt > /medpop/esp2/mesbah/datasets/ukbb/ukbb_v3/imp_v3/hapmap_v3.ukb_imp_chr1_v3


