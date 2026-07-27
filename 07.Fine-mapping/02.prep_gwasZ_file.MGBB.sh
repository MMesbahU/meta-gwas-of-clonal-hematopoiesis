#!/bin/bash

## while read CHR START END LEAD; do qsub -wd /broad/hptmp/mesbah/RefSeq/tmpdir -R y -l h_vmem=20G -l h_rt=02:00:00 -pe smp 1 -binding linear:1 -N prep_gwasZ.eaf01.mgbb53k /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/07.Fine-mapping/03.prep_gwasZ_file.sh /broad/hptmp/mesbah/dataset/mesbah_medpop_esp/GWAS_CHIP/MGBB53k/MultiANC/hasCH.chr1_22.MultiANC.regenie.tsv.gz /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/mgbb $CHR $START $END $LEAD 0.01 0.99; done < <(cut -f1-4 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/finemap_loci_500kb.hasCH.MultiANC.multianc_cohort6.txt)


# # mkdir -p /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/mgbb; 
# while read CHR START END LEAD; do bash /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/07.Fine-mapping/03.prep_gwasZ_file.sh /broad/hptmp/mesbah/dataset/mesbah_medpop_esp/GWAS_CHIP/MGBB53k/MultiANC/hasCH.chr1_22.MultiANC.regenie.tsv.gz /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/mgbb $CHR $START $END $LEAD; done < <(cut -f1-4 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/finemap_loci_500kb.hasCH.MultiANC.multianc_cohort6.txt) 1>>/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/mgbb/log_file.txt 2>>/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/mgbb/err_file.txt &
## running code
# while read CHR START END LEAD; do qsub -wd /broad/hptmp/mesbah/RefSeq/tmpdir -R y -l h_vmem=20G -l h_rt=02:00:00 -pe smp 1 -binding linear:1 -N prep_gwasZ.mgbb53k /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/07.Fine-mapping/03.prep_gwasZ_file.sh /broad/hptmp/mesbah/dataset/mesbah_medpop_esp/GWAS_CHIP/MGBB53k/MultiANC/hasCH.chr1_22.MultiANC.regenie.tsv.gz /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/mgbb $CHR $START $END $LEAD 0.01 0.99; done < <(cut -f1-4 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/finemap_loci_500kb.hasCH.MultiANC.multianc_cohort6.txt) 

# while read CHR START END LEAD; do qsub -wd /broad/hptmp/mesbah/RefSeq/tmpdir -R y -l h_vmem=20G -l h_rt=02:00:00 -pe smp 1 -binding linear:1 -N prep_gwasZ.eaf001.mgbb53k /medpop/esp2/mesbah/tools/meta-gwas-of-clonal-hematopoiesis/07.Fine-mapping/03.prep_gwasZ_file.sh /broad/hptmp/mesbah/dataset/mesbah_medpop_esp/GWAS_CHIP/MGBB53k/MultiANC/hasCH.chr1_22.MultiANC.regenie.tsv.gz /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/mgbb $CHR $START $END $LEAD 0.001 0.999; done < <(cut -f1-4 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/finemap_loci_500kb.hasCH.MultiANC.multianc_cohort6.txt)

# 1>>/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/mgbb/log_file.txt 2>>/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/Finemap/MMU/input/mgbb/err_file.txt &

#############
# /broad/hptmp/mesbah/dataset/mesbah_medpop_esp/GWAS_CHIP/MGBB53k/MultiANC/hasCH.chr1_22.MultiANC.regenie.tsv.gz | head -2
# SNPID	REF	ALT	AAF	BETA	SE	P	N	INFO
# chr1:10894:G:A	G	A	0.000978619	0.411001	0.573408	0.473516	50305	0.345706
# Z file The dataset.z file is a space-delimited text file and contains the GWAS summary statistics one SNP per line. It contains the mandatory column names in the following order.
# rsid column contains the SNP identifiers. The identifier can be a rsID number or a combination of chromosome name and genomic position (e.g. XXX:yyy)
# chromosome column contains the chromosome names. The chromosome names can be chosen freely with precomputed SNP correlations (e.g. 'X', '0X' or 'chrX')
# position column contains the base pair positions
# allele1 column contains the "first" allele of the SNPs. In SNPTEST this corresponds to 'allele_A', whereas BOLT-LMM uses 'ALLELE1'
# allele2 column contains the "second" allele of the SNPs. In SNPTEST this corresponds to 'allele_B', whereas BOLT-LMM uses 'ALLELE0'
# maf column contains the minor allele frequencies
# beta column contains the estimated effect sizes as given by GWAS software
# se column contains the standard errors of effect sizes as given by GWAS software
# flip optional column - see below 
# # SNPID	REF	ALT	AAF	BETA	SE	P	N	INFO
###############
GWAS=$1
OUTDir=$2
CHR=$3
START=$4
END=$5
LEAD=$6
EAF_min=${7}
EAF_max=${8}
OUT=${OUTDir}/${CHR}_${START}_${END}_${LEAD}.$(basename ${GWAS} ".regenie.tsv.gz").${EAF_min}_${EAF_max}.z
# rsid chromosome position allele1 allele2 maf beta se

zcat ${GWAS} | awk -v chr=$CHR -v start=$START -v end=$END -v eaf_min=$EAF_min -v eaf_max=$EAF_max '
BEGIN{
    OFS=" "
    print "rsid","chromosome","position","allele1","allele2","maf","beta","se"
    maxN=0
}

NR==1{
    for(i=1;i<=NF;i++){
        if($i=="SNPID") id=i
        if($i=="REF") ref=i
        if($i=="ALT") alt=i
        if($i=="AAF") maf=i
        if($i=="BETA") b=i
        if($i=="SE") se=i
        if($i=="P") pval=i
	if($i=="N") ns=i
    }
}

NR>1 {

    split($id,a,":")
    c=a[1]
    p=a[2]
    P=$pval+0
    MAF=$maf+0

    if(c==chr && p>=start && p<=end && P<0.5 && MAF>=eaf_min && MAF<=eaf_max){

        print $id,c,p,$alt,$ref,$maf,$b,$se

        # track max N in locus
        if(ns!=""){
            n=$ns+0
            if(n > maxN)
                maxN = n
        }
    }
}

END{
    print "## MAX_N_LOCI:", chr, start, end, maxN > "/dev/stderr"
}
' | awk 'BEGIN{OFS=" "} NR==1{print; next} {sub(/^chr/,"", $2); $6=($6>0.5)?1-$6:$6; print }' > ${OUT}

##

# SNP list
awk 'NR>1{print $1}' ${OUT} > $(dirname ${OUT})/rsid.$(basename ${OUT}).txt

