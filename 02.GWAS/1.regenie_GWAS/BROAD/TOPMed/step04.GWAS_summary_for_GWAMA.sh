#!/bin/bash


### AoU Workspace: "Detection of clonal hematopoiesis in 250k WGS [Dataset v7]"

# gsutil -m cp gs://fc-secure-a1a19750-8c62-44db-bc3e-d2baed266066/Meta4/sum4_cohort/STRATIFIED/topmed/*.gz /home/jupyter/metagwas/topmed/

## Prepare summary for gwama:
# Note: Alt="Genome Reference Allele"
## Effect allele is the true REF allele
## Name    Chr     Pos     Ref     Alt     Trait   Cohort  Model   Effect  LCI_Effect      UCI_Effect      Pval    AAF     Num_Cases  Cases_Ref       Cases_Het       Cases_Alt       Num_Controls    Controls_Ref    Controls_Het    Controls_Alt    Info
# rs886376966     1       133857  A       G       hasCHvaf02      AMR.TOPMed74k   ADD-WGR-FIRTH   0.414975        0.0581053     2.96366 0.380567        0.998445        398     0       2       396     7639    0       23      7616    REGENIE_BETA=-0.879536;REGENIE_SE=1.003060;INFO=1.000000;MAC=25.000000;SCORE=0.874178;SKATV=0.993908;LOG10P=0.419569
while read ANC; do mv ~/metagwas/topmed/${ANC}.TOPMed74k.chr1_22.hasCHvaf02.regenie.gz ~/metagwas/topmed/${ANC}.TOPMed74k.chr1_22.hasCH.regenie.gz; done < <(echo -e "AMR\nAFR\nEUR\nFemale\nMale\nMultiANC")

while read ANC; do while read phenos; do echo -e "SNPID\tRSID\tREF\tALT\tAAF\tBETA\tSE\tP\tN" | gzip -c > ~/metagwas/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${phenos}.regenie.tsv.gz && zcat ~/metagwas/topmed/${ANC}.TOPMed74k.chr1_22.has${phenos}.regenie.gz | sed '1d' | sed -e 's:;:\t:g' -e 's:REGENIE_BETA=::g' -e 's:REGENIE_SE=::g' | awk '{print "chr"$2":"$3":"$5":"$4"\t"$1"\t"$4"\t"$5"\t"$13"\t"$22"\t"$23"\t"$12"\t"($14+$18)}' | gzip -c >> ~/metagwas/topmed/GWAMA.${ANC}.TOPMed74k.chr1_22.has${phenos}.regenie.tsv.gz; done< <(echo -e "CH\nDNMT3A\nTET2"); done < <(echo -e "AFR\nAMR\nEUR\nFemale\nMale\nMultiANC") &


