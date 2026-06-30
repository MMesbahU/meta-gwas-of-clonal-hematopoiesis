
#!/bin/sh
#$ -wd /medpop/esp2/btruong/Projects/logjobs
#$ -j y
#$ -pe smp 8 -R y -binding linear:8
#$ -l h_vmem=8G
#$ -l h_rt=20:00:00
#$ -N magma


source /broad/software/scripts/useuse
source ~/.my.bashrc
use GCC-5.2

reuse Python-3.6
reuse Anaconda
reuse Anaconda3


cd ${wdir}


/medpop/esp2/btruong/Tools/magma --annotate --snp-loc ${trait}.snploc --gene-loc ${geneloc} --out magma_annot_${trait}

/medpop/esp2/btruong/Tools/magma \
  --bfile /broad/hptmp/btruong/HDP/pops/g1000_eur \
  --gene-annot magma_annot_${trait}.genes.annot \
  --pval ${trait}.snploc ncol=N \
  --gene-model snp-wise=mean \
  --out magma_${trait}_1

