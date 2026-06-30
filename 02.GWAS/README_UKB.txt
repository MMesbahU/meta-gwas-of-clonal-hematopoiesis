# UKB analysis

1. CHIP GWAS: (separate GWAS for 200k and 250k  WES samples) 
	1.1 NULL step in GCP 
	1.2 SVAS in Broad HPC
2. Prepare GWAS Summary Stats using "/medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/prepGWASsummary.UKB.sh"

3. LiftOver UKB hg37 loci to hg38 for meta-analysis '/medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/01.PrepSummary/02.prepUKBB_summary.sh'

4. Perform Meta-Analysis using GWAMA and add 'N' from METAL summary stats

