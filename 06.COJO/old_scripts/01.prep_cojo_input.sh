#!/bin/bash

## for files in $(ls -l /medpop/esp/mesbah/GWAS_CHIP/meta_gwas/GWAMA_VCF/hg37/lifted_hg37.GWAMA.m*.vcf.gz| awk '{print $NF}'); do qsub -R y -wd /broad/hptmp/mesbah/ukb_chip/meta_gwas/cojo/tmpdir -N prep_cojo -l h_rt=20:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/04.COJO/01.prep_cojo_input.sh ${files} /broad/hptmp/mesbah/ukb_chip/meta_gwas/cojo/cojo_input; done

#### 
source /broad/software/scripts/useuse

use Bcftools
use Tabix


### Read command line inputs
gwas_vcf=${1}

cojo_input_dir=${2}

### Prepare COJO input from VCF file
cojo_input=${cojo_input_dir}/summary_for_cojo.$(basename ${gwas_vcf} ".vcf.gz").tsv

echo -e "SNP\tA1\tA2\tfreq\tb\tse\tp\tN" > ${cojo_input} 

bcftools view  -i 'EAF>=0.01 & EAF<=0.99' ${gwas_vcf} | bcftools query -f '%ID\t%ALT\t%REF\t%INFO/EAF\t%INFO/BETA\t%INFO/SE\t%INFO/pvalue\t%INFO/N_METAL\n' >> ${cojo_input}

