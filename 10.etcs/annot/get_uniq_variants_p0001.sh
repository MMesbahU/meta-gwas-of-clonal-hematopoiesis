while read files; do zcat ${files} | awk 'NR>1{print $1}' | sort -V | uniq >> /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/summary/any_varID_hg38_p0001.list; done < <(ls /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/summary/pval_0.0001.GWAMA.chr1_22.has*.tsv.gz | awk '{print $NF}') &

## 
cat /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/summary/any_varID_hg38_p0001.list | sort -V | uniq > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/summary/unique_varID_hg38_p0001.list

gzip /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/summary/*.list

## Prep VCF
inVCF=/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/summary/unique_varID_hg38_p0001.vcf
echo -e '##fileformat=VCFv4.2' > ${inVCF}
echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' >> ${inVCF}

zcat /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/summary/unique_varID_hg38_p0001.list.gz | sed 's:\::\t:g' | awk '{print $1"\t"$2"\t"$1":"$2":"$3":"$4"\t"$3"\t"$4"\t.\t.\t."}' >> ${inVCF}

gzip ${inVCF}

### Annovar
cp /medpop/esp2/projects/software/annovar/humandb/hg38_* /broad/hptmp/mesbah/dataset/annovar/humandb/ &

gzip -d /broad/hptmp/mesbah/dataset/annovar/humandb/hg38_*.gz


