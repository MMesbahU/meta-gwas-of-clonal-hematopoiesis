### GWAS significant SNPs, P<=1.67e-8
while read lines; do zcat ${lines} | head -1 > /medpop/esp2/mesbah/projects/Meta_GWAS/sumstats/varSig/p167e8.$(basename ${lines} ".tsv.gz").tsv && zcat ${lines} | awk '$12<=1.67e-8' >> /medpop/esp2/mesbah/projects/Meta_GWAS/sumstats/varSig/p167e8.$(basename ${lines} ".tsv.gz").tsv ;done < <(ls -lh /medpop/esp2/mesbah/projects/Meta_GWAS/sumstats/lifted_hg37.GWAMA.meta_*.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.tsv.gz | awk '{print $NF}') &


## zgrep -wE 'chr1:226422811:C:G|chr3:160368267:G:T|chr4:104801703:A:G|chr4:104807825:G:A|chr5:1285859:C:A|chr6:109360679:T:C|chr11:108490851:G:A|chr12:26437841:G:A|chr17:57372795:T:G|chr18:44569779:G:A|chr22:28695868:AG:A' /medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/topmed2019.all_samples.hasCHIP.tsv.gz
## CHIP
while read files; do zcat ${files} | head -1 > /medpop/esp2/mesbah/projects/Meta_GWAS/sumstats/leadVars/$(basename ${files} ".tsv.gz").tsv && zgrep -wE 'chr1:226422811:C:G|chr3:160368267:G:T|chr4:104801703:A:G|chr4:104807825:G:A|chr5:1285859:C:A|chr6:109360679:T:C|chr11:108490851:G:A|chr12:26437841:G:A|chr17:57372795:T:G|chr18:44569779:G:A|chr22:28695868:AG:A' ${files} >> /medpop/esp2/mesbah/projects/Meta_GWAS/sumstats/leadVars/$(basename ${files} ".tsv.gz").tsv; done </medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_GWAMA/Input_GWAMA.allsamples.topmed_ukb_mgb_biovu.CHIP_GWAS.txt

## DNMT3A
while read files; do zcat ${files} | head -1 > /medpop/esp2/mesbah/projects/Meta_GWAS/sumstats/leadVars/$(basename ${files} ".tsv.gz").tsv && zgrep -wE 'chr1:226422811:C:G|chr3:160368267:G:T|chr3:160364283:T:C|chr4:104801703:A:G|chr4:104807825:G:A|chr5:1285859:C:A|chr6:109478720:C:T|chr6:109360679:T:C|chr11:108484147:C:G|chr11:108490851:G:A|chr12:26422298:GAAT:G|chr12:26437841:G:A|chr14:95714358:G:T|chr18:44581678:G:T|chr18:44569779:G:A|chr22:28695868:AG:A' ${files} >> /medpop/esp2/mesbah/projects/Meta_GWAS/sumstats/leadVars/$(basename ${files} ".tsv.gz").tsv; done </medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_GWAMA/Input_GWAMA.allsamples.topmed_ukb_mgb_biovu.DNMT3A_GWAS.txt

## TET2
while read files; do zcat ${files} | head -1 > /medpop/esp2/mesbah/projects/Meta_GWAS/sumstats/leadVars/$(basename ${files} ".tsv.gz").tsv && zgrep -wE 'chr1:11321655:C:G|chr5:1286401:C:A|chr5:1285859:C:A|chr6:131782542:T:C|chr7:42286311:G:A|chr14:95727218:C:T|chr14:95714358:G:T|chr16:87972320:A:T' ${files} >> /medpop/esp2/mesbah/projects/Meta_GWAS/sumstats/leadVars/$(basename ${files} ".tsv.gz").tsv; done </medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_GWAMA/Input_GWAMA.allsamples.topmed_ukb_mgb_biovu.TET2_GWAS.txt


