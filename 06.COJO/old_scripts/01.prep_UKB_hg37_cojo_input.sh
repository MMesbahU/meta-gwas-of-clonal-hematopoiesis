#!/bin/bash

## for files in $(ls -l /medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/chr1_22.has*.21Aug2021_ukbEUR.tsv.gz | awk '{print $NF}'); do qsub -R y -wd /broad/hptmp/mesbah/ukb_chip/meta_gwas/cojo/tmpdir -N prep_cojo -l h_rt=20:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.COJO/01.prep_cojo_input.sh ${files} /broad/hptmp/mesbah/ukb_chip/meta_gwas/cojo/cojo_input; done

#### 
# source /broad/software/scripts/useuse

# use Bcftools
# use Tabix


### Read command line inputs
gwas_tsv=${1}

cojo_input_dir=${2}

### Prepare COJO input from VCF file
cojo_input=${cojo_input_dir}/summary_for_cojo.$(basename ${gwas_tsv} ".gz")

echo -e "SNP\tA1\tA2\tfreq\tb\tse\tp\tN" > ${cojo_input} 

zcat ${gwas_tsv} | awk 'NR>1 && $13>=0.01 && $13<=0.99 && $24>=0.9{print $1"\t"$5"\t"$4"\t"$13"\t"$22"\t"$23"\t"$12"\t"($14+$18)}' >> ${cojo_input}

