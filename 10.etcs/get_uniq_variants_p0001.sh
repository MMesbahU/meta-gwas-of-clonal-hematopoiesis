while read files; do zcat ${files} | awk 'NR>1{print $1}' | sort -V | uniq >> /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/summary/any_varID_hg38_p0001.list; done < <(ls /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/summary/pval_0.0001.GWAMA.chr1_22.has*.tsv.gz | awk '{print $NF}') &

## 
cat /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/summary/any_varID_hg38_p0001.list | sort -V | uniq > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/summary/unique_varID_hg38_p0001.list

## Prep VCF

