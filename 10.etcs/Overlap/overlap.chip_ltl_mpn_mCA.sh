### mCA GWAS in UKB: Zekavat, S.M. et al. Hematopoietic mosaic chromosomal alterations increase the risk for diverse types of infection. Nat Med 27, 1012-1024 (2021)
zcat /medpop/esp/mesbah/GWAS_CHIP/mCA/Maryam.expanded_mCA_summary_stats.N444199.gz | awk '(NR==1){print $0}(NR>1 && $6>=7.30103){print $0}' > /medpop/esp/mesbah/GWAS_CHIP/meta_gwas/Overlaps/p5e8.Maryam.expanded_mCA_summary_stats.N444199.tsv &

## LTL gwas in UKB: Codd, V. et al. Polygenic basis and biomedical consequences of telomere length variation. Nat Genet 53, 1425-1433 (2021).  
zcat /medpop/esp/mesbah/GWAS_CHIP/ltl_gwas/UKB_telomere_gwas_summarystats.tsv.gz | awk '(NR==1){print $0}(NR>1 && $2<=5e-8){print $0}' > /medpop/esp/mesbah/GWAS_CHIP/meta_gwas/Overlaps/p5e8.UKB_telomere_gwas_summarystats.tsv &

## MPN: Bao, E.L. et al. Inherited myeloproliferative neoplasm risk affects haematopoietic stem cells. Nature 586, 769-775 (2020). 
zcat /medpop/esp/mesbah/GWAS_CHIP/MPN_summary/MPN_metaGWAS_sumstats.tsv.gz | awk '(NR==1){print $0}(NR>1 && $9<=5e-8){print $0}' > /medpop/esp/mesbah/GWAS_CHIP/meta_gwas/Overlaps/p5e8.MPN_metaGWAS_sumstats.tsv &

## CHIP 
zcat /medpop/esp2/mesbah/projects/Meta_GWAS/sumstats/chr1_22.lifted_hg37.GWAMA.meta_CHIP.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.tsv.gz | awk '(NR==1){print $0}(NR>1 && $12<=5e-8){print $0}' > /medpop/esp/mesbah/GWAS_CHIP/meta_gwas/Overlaps/p5e8.chr1_22.lifted_hg37.GWAMA.meta_CHIP.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.tsv &

## DNMT3A
zcat /medpop/esp2/mesbah/projects/Meta_GWAS/sumstats/chr1_22.lifted_hg37.GWAMA.meta_DNMT3A.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.tsv.gz | awk '(NR==1){print $0}(NR>1 && $12<=5e-8){print $0}' > /medpop/esp/mesbah/GWAS_CHIP/meta_gwas/Overlaps/p5e8.chr1_22.lifted_hg37.GWAMA.meta_DNMT3A.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.tsv &

## TET2 
zcat /medpop/esp2/mesbah/projects/Meta_GWAS/sumstats/chr1_22.lifted_hg37.GWAMA.meta_TET2.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.tsv.gz | awk '(NR==1){print $0}(NR>1 && $12<=5e-8){print $0}' > /medpop/esp/mesbah/GWAS_CHIP/meta_gwas/Overlaps/p5e8.chr1_22.lifted_hg37.GWAMA.meta_TET2.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.tsv &

## Multi-ancestry
intersectBed -a <( awk 'NR>1{print $2"\t"$3-1"\t"$3"\t"$0}' /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/kar.chip.GCST90102618_buildGRCh37.p5e8.tsv ) -b <( awk 'NR>1{print $1"\t"$2-1"\t"$2"\t"$0}' /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/summary_p5e8.CHIP.tsv) -wao | awk '$NF>0' > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/bedintersect.chip.kar.tsv

intersectBed -a <( awk 'NR>1{print $2"\t"$3-1"\t"$3"\t"$0}' /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/kar.dnmt3a.GCST90102619.p5e8.tsv ) -b <( awk 'NR>1{print $1"\t"$2-1"\t"$2"\t"$0}' /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/summary_p5e8.DNMT3A.tsv) -wao | awk '$NF>0' > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/bedintersect.dnmt3a.kar.tsv

intersectBed -a <( awk 'NR>1{print $2"\t"$3-1"\t"$3"\t"$0}' /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/kar.tet2.GCST90102620.p5e8.tsv ) -b <( awk 'NR>1{print $1"\t"$2-1"\t"$2"\t"$0}' /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/summary_p5e8.TET2.tsv) -wao | awk '$NF>0' > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/bedintersect.tet2.kar.tsv

## EUR only 
intersectBed -a <( awk 'NR>1{print $2"\t"$3-1"\t"$3"\t"$0}' /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/kar.chip.GCST90102618_buildGRCh37.p5e8.tsv ) -b <( awk 'NR>1{print $1"\t"$2-1"\t"$2"\t"$0}' /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/eur_summary_p5e8.CHIP.tsv) -wao | awk '$NF>0' > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/bedintersect.eur_chip.kar.tsv

intersectBed -a <( awk 'NR>1{print $2"\t"$3-1"\t"$3"\t"$0}' /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/kar.dnmt3a.GCST90102619.p5e8.tsv ) -b <( awk 'NR>1{print $1"\t"$2-1"\t"$2"\t"$0}' /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/eur_summary_p5e8.DNMT3A.tsv) -wao | awk '$NF>0' > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/bedintersect.eur_dnmt3a.kar.tsv

intersectBed -a <( awk 'NR>1{print $2"\t"$3-1"\t"$3"\t"$0}' /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/kar.tet2.GCST90102620.p5e8.tsv ) -b <( awk 'NR>1{print $1"\t"$2-1"\t"$2"\t"$0}' /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/eur_summary_p5e8.TET2.tsv) -wao | awk '$NF>0' > /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/bedintersect.eur_tet2.kar.tsv

