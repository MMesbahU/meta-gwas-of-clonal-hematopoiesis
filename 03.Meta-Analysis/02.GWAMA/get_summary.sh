head -1 /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.hasDDR.AFR.ukbb450k_AoU250k_TOPMed72k_MGBB53k_ppm1dMCPS136k.out > /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/out/p5e8.hasDDR.AFR_MCPS136k.tsv; awk '(NR>1 && $10<5e-8){print $0}' /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.chr1_22.hasDDR.AFR.ukbb450k_AoU250k_TOPMed72k_MGBB53k_ppm1dMCPS136k.out | sort -V -k1 >> /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/out/p5e8.hasDDR.AFR_MCPS136k.tsv

## 
mkdir -p /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/out/boost_mcps

while read lines; do zcat ${lines} | head -1 > /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/out/boost_mcps/$(basename ${lines} ".out.gz").tsv; zcat ${lines} | awk '(NR>1 && $10<5e-8){print $0}' | sort -V -k1 >> /broad/hptmp/mesbah/dataset/ch_gwas/GWAMA/out/boost_mcps/$(basename ${lines} ".out.gz").tsv; done < <(ls /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/GWAMA/out/GWAMA.*CPS136k.out.gz | awk '{print $NF}') &





