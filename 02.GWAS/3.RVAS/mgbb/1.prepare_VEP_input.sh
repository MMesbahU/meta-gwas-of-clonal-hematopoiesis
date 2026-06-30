#!/bin/bash

## for chr in {1..22}; do bash /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/02_2.RVAS/prepare_VEP_input.sh /medpop/esp2/projects/MGB_Biobank/exome/53K_WES/release/annotation/MGB_53K_WES_genotype_variant_sample_QCed_chr${chr}.vep.gz /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/rvas/mask_files/chr${chr}.MGB_53K_WES.vep.hc_LOF.annotation; done
in_vep=${1}
out_regenie_annot=${2}
zgrep -v '^##' ${in_vep}  | \
awk '
BEGIN{FS=OFS="\t"}
NR==1{
  for(i=1;i<=NF;i++){
  	if($i=="#Uploaded_variation") vid=i
        if($i=="Gene") gene=i
	if($i=="IMPACT") impact=i
	if($i=="LoF") lof=i
     }
     next   
}
$impact=="HIGH" && $lof=="HC" {
	print $vid,$gene,"LoF"
}'  > ${out_regenie_annot}

