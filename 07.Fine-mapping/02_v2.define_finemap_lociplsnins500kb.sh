#!/bin/bash

# bash /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/07.Fine-mapping/02_v2.define_finemap_lociplsnins500kb.sh /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/lead_snp.significant.hasCH.MultiANC.multianc_cohort6.txt /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/finemap_loci_500kb.hasCH.MultiANC.multianc_cohort6.txt


#!/bin/bash

INPUT=$1
OUTPUT=$2

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
