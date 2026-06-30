#!/bin/bash

## qsub -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas/tmpdir -R y -l h_vmem=20G -l h_rt=20:00:00 -pe smp 1 -binding linear:1 -N topmed_gwas_summary /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/BROAD/TOPMed/step04.GWAS_summary_for_GWAMA.local.sh /broad/hptmp/mesbah/dataset/ch_gwas/topmed/gwas/2024 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas

### AoU Workspace: "Detection of clonal hematopoiesis in 250k WGS [Dataset v7]"

# gsutil -m cp gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/Meta4/sum4_cohort/STRATIFIED/topmed/*.gz /home/jupyter/metagwas/topmed/

## Prepare summary for gwama:
# Note: Alt="Genome Reference Allele"
## Effect allele is the true REF allele
## Name    Chr     Pos     Ref     Alt     Trait   Cohort  Model   Effect  LCI_Effect      UCI_Effect      Pval    AAF     Num_Cases  Cases_Ref       Cases_Het       Cases_Alt       Num_Controls    Controls_Ref    Controls_Het    Controls_Alt    Info
# rs886376966     1       133857  A       G       hasCHvaf02      AMR.TOPMed74k   ADD-WGR-FIRTH   0.414975        0.0581053     2.96366 0.380567        0.998445        398     0       2       396     7639    0       23      7616    REGENIE_BETA=-0.879536;REGENIE_SE=1.003060;INFO=1.000000;MAC=25.000000;SCORE=0.874178;SKATV=0.993908;LOG10P=0.419569
##
outdir=${1} # /broad/hptmp/mesbah/dataset/ch_gwas/topmed/gwas/2024

indir=${2} # /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/topmed_v10/gwas
##
while read ANC
do 
	mv ${indir}/${ANC}/${ANC}.TOPMed74k.chr1_22.hasCHvaf02.regenie.gz ${indir}/${ANC}/${ANC}.TOPMed74k.chr1_22.hasCH.regenie.gz
done < <(echo -e "AMR\nAFR\nEUR\nFemale\nMale\nMultiANC")

##
while read ANC
do 
	while read phenos
	do 
		
		rm ${outdir}/GWAMA.${ANC}.TOPMed74k.chr1_22.has${phenos}.regenie.tsv.gz

		echo -e "SNPID\tRSID\tREF\tALT\tAAF\tBETA\tSE\tP\tN" > ${outdir}/GWAMA.${ANC}.TOPMed74k.chr1_22.has${phenos}.regenie.tsv 
			
		zcat ${indir}/${ANC}/${ANC}.TOPMed74k.chr1_22.has${phenos}.regenie.gz | sed '1d' | sed -e 's:;:\t:g' -e 's:REGENIE_BETA=::g' -e 's:REGENIE_SE=::g' | awk '{print "chr"$2":"$3":"$5":"$4"\t"$1"\t"$4"\t"$5"\t"$13"\t"$22"\t"$23"\t"$12"\t"($14+$18)}' >> ${outdir}/GWAMA.${ANC}.TOPMed74k.chr1_22.has${phenos}.regenie.tsv 
		
		gzip -f ${outdir}/GWAMA.${ANC}.TOPMed74k.chr1_22.has${phenos}.regenie.tsv 

	done< <(echo -e "CH\nCHvaf10\nDNMT3A\nTET2\nASXL1\nSF\nDDR")
	
done < <(echo -e "AFR\nAMR\nEUR\nFemale\nMale\nMultiANC")



