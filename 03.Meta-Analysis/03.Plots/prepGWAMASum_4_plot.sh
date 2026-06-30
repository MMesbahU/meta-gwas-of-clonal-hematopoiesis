#!/bin/bash

# To Run: bash script.sh in_summary.gz out_summary.tsv
## while read pheno; do qsub -R y -wd /broad/hptmp/mesbah/dataset/ch_gwas/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=20G -l h_rt=10:00:00 -N prep.${pheno}_sum_for_plot /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/prepGWAMASum_4_plot.sh /broad/hptmp/mesbah/dataset/ch_gwas/metaGWAS.has${pheno}.topmed74k_topImpukb200k250k_aou98k_mgbb53k_biovu54k.out.gz /broad/hptmp/mesbah/dataset/ch_gwas/for_plot.metaGWAS.has${pheno}.topmed74k_topImpukb200k250k_aou98k_mgbb53k_biovu54k.nStd3pls.tsv; done < <(echo -e "CH\nDNMT3A\nTET2")

# GWAMA header
# rs_number	reference_allele	other_allele	eaf	beta	se	beta_95L	beta_95U	z	p-value	_-log10_p-value	q_statistic	q_p-value	i2	n_studies	n_samples	effects
# chr1:19776:A:G	G	A	0.001887	-0.378736	0.372975	-1.109767	0.352296	-1.015445	0.309880	0.508807	-0.000000	1.000000	0.000000	1	74974	-?????

# manhattan plot file with 4 columns
# c("SNP","Chr", "POS", "P")

gwas_sum=${1}

plot_sum=${2}

echo -e "SNP\tChr\tPOS\tP" > ${plot_sum}

zcat ${gwas_sum} | awk 'NR>1 && $4>=0.001 && $4<=0.999 && $15>2 {print $1"\t"$10}' | sed 's:\::\t:g' | awk '{print $1":"$2":"$3":"$4"\t"$1"\t"$2"\t"$5}' >> ${plot_sum}

gzip -f ${plot_sum}


