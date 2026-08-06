#!/bin/bash

# bash /home/jupyter/meta-gwas-of-clonal-hematopoiesis/02.GWAS/3.RVAS/aou/1.prepare_VEP_input.sh /home/jupyter/exome_v7/annotation/all/AoUv7_exome.genotypeQCed.variantQCed.sampleQCed.vep.gz /home/jupyter/exome_v7/mask_files/chr1_22.AoUv7_exome.genotypeQCed.variantQCed.sampleQCed.vep.hc_LOF.annotation

## for chr in {1..22}; do bash /home/jupyter/meta-gwas-of-clonal-hematopoiesis/02.GWAS/3.RVAS/aou/1.prepare_VEP_input.sh /home/jupyter/exome_v7/annotation/all/AoUv7_exome.genotypeQCed.variantQCed.sampleQCed.vep.gz /home/jupyter/exome_v7/mask_files/chr$chr.AoUv7_exome.genotypeQCed.variantQCed.sampleQCed.vep.hc_LOF.annotation $chr; done

#
# #CHR	POS	ID	REF	ALT	Uploaded_variation	Location	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	IMPACT	DISTANCE	STRAND	FLAGS	SYMBOL	SYMBOL_SOURCE	HGNC_ID	CANONICAL	CLIN_SIG	SOMATIC	PHENO	LoF	LoF_filter	LoF_flags	LoF_info

#

in_vep=${1}

out_regenie_annot=${2}

# chr=${3}
# 
zcat ${in_vep}  | \
awk '
BEGIN{FS=OFS="\t"}
NR==1{
  for(i=1;i<=NF;i++){
	if($i=="ID") vid=i
        if($i=="Gene") gene=i
	if($i=="SYMBOL") gene_SYMBOL=i
	if($i=="IMPACT") impact=i
	if($i=="LoF") lof=i
     }
     next   
}
$impact=="HIGH" && $lof=="HC" {
	print $vid,$gene":"$gene_SYMBOL,"LoF"
}' | awk '!seen[$1]++' > ${out_regenie_annot}

## for chr in {1..22} X; do  awk -v chrom="chr${chr}:" '$1 ~ chrom' /home/jupyter/exome_v7/mask_files/chr1_22.AoUv7_exome.genotypeQCed.variantQCed.sampleQCed.vep.hc_LOF.annotation > /home/jupyter/exome_v7/mask_files/chr${chr}.AoUv7_exome.genotypeQCed.variantQCed.sampleQCed.vep.hc_LOF.annotation; done


