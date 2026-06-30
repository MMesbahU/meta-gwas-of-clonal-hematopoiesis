#!/bin/bash


# 10:100344209:C:T
tmpVCF=/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/mash/2025_run/varID.vcf

echo -e '##fileformat=VCFv4.2' > ${tmpVCF}

echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' >> ${tmpVCF}

zcat /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/mash/2025_run/chip_dnmt_tet_ltl_mCA_mpn.lfsr05.oct2024.csv.gz | awk -F "," 'NR>1{print $1}' | sort -V | uniq | sed 's:\::\t:g' | awk '{print $1"\t"$2"\t"$1":"$2":"$3":"$4"\t"$3"\t"$4"\t.\t.\t."}' >> ${tmpVCF}

## Annotate
gzip -d /broad/hptmp/mesbah/dataset/annovar/humandb/hg19_*.gz

#

##
