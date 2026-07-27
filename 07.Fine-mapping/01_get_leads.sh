#!/bin/bash

## Multi-ancestry:
#afreq=0.01; while read files; do bash /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/07.Fine-mapping/01_get_leads.sh ${files} /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/lead_snp_p5e8_eaf${afreq}.$(basename ${files} ".out.gz").tsv ${afreq}; done < <(ls /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has*.MultiANC.ukbb200k_ukb250k*out.gz | awk '{print $NF}') &

## AFR:
# afreq=0.01; while read files; do bash /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/07.Fine-mapping/01_get_leads.sh ${files} /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/lead_snp_p5e8_eaf${afreq}.$(basename ${files} ".out.gz").tsv ${afreq}; done < <(ls /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has*.AFR.ukbb450k*.out.gz | awk '{print $NF}') &

# AMR:
# afreq=0.01; while read files; do bash /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/07.Fine-mapping/01_get_leads.sh ${files} /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/lead_snp_p5e8_eaf${afreq}.$(basename ${files} ".out.gz").tsv ${afreq}; done < <(ls /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has*.AMR.ukbb450k*.out.gz | awk '{print $NF}') &
# EUR:
# afreq=0.01; while read files; do bash /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/07.Fine-mapping/01_get_leads.sh ${files} /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/lead_snp_p5e8_eaf${afreq}.$(basename ${files} ".out.gz").tsv ${afreq}; done < <(ls /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has*.EUR.ukbb450k*.out.gz | awk '{print $NF}') &

# Male:
# afreq=0.01; while read files; do bash /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/07.Fine-mapping/01_get_leads.sh ${files} /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/lead_snp_p5e8_eaf${afreq}.$(basename ${files} ".out.gz").tsv ${afreq}; done < <(ls /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has*.Male.ukbb450k*.out.gz | awk '{print $NF}') &

# Female:
# afreq=0.01; while read files; do bash /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/07.Fine-mapping/01_get_leads.sh ${files} /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/lead_snp_p5e8_eaf${afreq}.$(basename ${files} ".out.gz").tsv ${afreq}; done < <(ls /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has*.Female.ukbb450k*.out.gz | awk '{print $NF}') &

#### test run ####
## bash /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/07.Fine-mapping/01_get_leads.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.hasCH.MultiANC.ukbb200k_ukb250k_AoU250k_TOPMed72k_MGBB53k_BioVU54k_MCPS136k.out.gz /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/lead_snp.significant.txt


GWAS=$1
OUT=$2
allele_freq=$3
##
zcat ${GWAS} \
| awk  -v afeq=${allele_freq} -v OFS="\t" '
NR==1{
    for(i=1;i<=NF;i++){
        if($i=="rs_number") id=i
        if($i=="eaf") eaf=i
        if($i=="p-value") pval=i
        if($i=="n_studies") nstudies=i
	if($i=="n_samples") nsample=i
      }
}
NR>1{
    split($id, a, ":")
    chr=a[1]
    pos=a[2]
    variant=$id

    P=$pval+0
    EAF=$eaf+0
    N=$nsample+0
    n_st=$nstudies+0

    if(P <= 5e-8 && EAF > 0.001 && EAF <0.999 && n_st >=2) {
        print chr, pos, variant, P, EAF, N
	}
}' \
| sort -k3 -V \
> ${OUT}

