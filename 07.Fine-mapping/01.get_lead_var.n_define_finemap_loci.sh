#!/bin/bash
## Multi-ancestry:
# afreq=0.01; while read files; do qsub -wd /broad/hptmp/mesbah/RefSeq/tmpdir -R y -l h_vmem=20G -l h_rt=02:00:00 -pe smp 1 -binding linear:1 -N prep_gwas_loci.multianc /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/07.Fine-mapping/01.get_lead_var.n_define_finemap_loci.sh ${files} /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/lead_snp_p5e8_eaf${afreq}.$(basename ${files} ".out.gz").tsv ${afreq}; done < <(ls /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has*.MultiANC.ukbb200k_ukb250k*out.gz | awk '{print $NF}')

## AFR:
# afreq=0.01; while read files; do qsub -wd /broad/hptmp/mesbah/RefSeq/tmpdir -R y -l h_vmem=20G -l h_rt=02:00:00 -pe smp 1 -binding linear:1 -N prep_gwas_loci.afr /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/07.Fine-mapping/01.get_lead_var.n_define_finemap_loci.sh ${files} /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/lead_snp_p5e8_eaf${afreq}.$(basename ${files} ".out.gz").tsv ${afreq}; done < <(ls /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has*.AFR.ukbb450k*.out.gz | awk '{print $NF}')

# AMR:
# afreq=0.01; while read files; do qsub -wd /broad/hptmp/mesbah/RefSeq/tmpdir -R y -l h_vmem=20G -l h_rt=02:00:00 -pe smp 1 -binding linear:1 -N prep_gwas_loci.amr /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/07.Fine-mapping/01.get_lead_var.n_define_finemap_loci.sh ${files} /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/lead_snp_p5e8_eaf${afreq}.$(basename ${files} ".out.gz").tsv ${afreq}; done < <(ls /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has*.AMR.ukbb450k*.out.gz | awk '{print $NF}') 

# EUR:
# afreq=0.01; while read files; do qsub -wd /broad/hptmp/mesbah/RefSeq/tmpdir -R y -l h_vmem=20G -l h_rt=02:00:00 -pe smp 1 -binding linear:1 -N prep_gwas_loci.eur /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/07.Fine-mapping/01.get_lead_var.n_define_finemap_loci.sh ${files} /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/lead_snp_p5e8_eaf${afreq}.$(basename ${files} ".out.gz").tsv ${afreq}; done < <(ls /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has*.EUR.ukbb450k*.out.gz | awk '{print $NF}')

# Male:
# afreq=0.01; while read files; do qsub -wd /broad/hptmp/mesbah/RefSeq/tmpdir -R y -l h_vmem=20G -l h_rt=02:00:00 -pe smp 1 -binding linear:1 -N prep_gwas_loci.male /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/07.Fine-mapping/01.get_lead_var.n_define_finemap_loci.sh ${files} /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/lead_snp_p5e8_eaf${afreq}.$(basename ${files} ".out.gz").tsv ${afreq}; done < <(ls /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has*.Male.ukbb450k*.out.gz | awk '{print $NF}')

# Female:
# afreq=0.01; while read files; do qsub -wd /broad/hptmp/mesbah/RefSeq/tmpdir -R y -l h_vmem=20G -l h_rt=02:00:00 -pe smp 1 -binding linear:1 -N prep_gwas_loci.female /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/07.Fine-mapping/01.get_lead_var.n_define_finemap_loci.sh ${files} /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/lead_snp_p5e8_eaf${afreq}.$(basename ${files} ".out.gz").tsv ${afreq}; done < <(ls /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.has*.Female.ukbb450k*.out.gz | awk '{print $NF}')

#### test run ####
## bash /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/07.Fine-mapping/01_get_leads.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.hasCH.MultiANC.ukbb200k_ukb250k_AoU250k_TOPMed72k_MGBB53k_BioVU54k_MCPS136k.out.gz /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/lead_snp.significant.txt

#######
GWAS=$1
OUT=$2
allele_freq=$3

# inputs for 2nd task
INPUT=${OUT}
OUTPUT=$(dirname ${OUT})/finemap_loci_500kb.$(basename ${OUT})

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

    if(P <= 5e-8 && EAF > afeq && EAF < (1 - afeq) && n_st >=2) {
        print chr, pos, variant, P, EAF, N
	}
}' \
| sort -k3 -V \
> ${OUT}

### Lead variants and define  Loci ###
# INPUT=$1
# OUTPUT=$2

awk '
BEGIN{
    OFS="\t"
}

NR==1{
    chr=$1
    start=$2
    end=$2

    min_p=1
    lead=$3

    n=1
    next
}

{
    c=$1
    pos=$2
    snp=$3

    # new locus condition
    if(c!=chr || (pos-end)>500000){

        # finalize previous locus
        fm_start=start-500000
        fm_end=end+500000

        if(fm_start<1) fm_start=1

        print chr,
              fm_start,
              fm_end,
              lead,
              n,
              start,
              end

        # reset locus
        chr=c
        start=pos
        end=pos
        min_p=1
        lead=snp
        n=1

    } else {

        # extend locus
        if(pos > end) end=pos

        n++
    }

    # update lead SNP by minimum p-value (ASSUMES 4th column = P-value)
    p=$4+0

    if(p < min_p){
        min_p=p
        lead=snp
    }
}

END{
    fm_start=start-500000
    fm_end=end+500000

    if(fm_start<1) fm_start=1

    print chr,
          fm_start,
          fm_end,
          lead,
          n,
          start,
          end
}
' ${INPUT} > ${OUTPUT}


