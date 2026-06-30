## R v 4.1.1 (2021-08-10)

library(data.table)

## MPN GWAS: source 'Bao, E. L. et al. Inherited myeloproliferative neoplasm risk affects haematopoietic stem cells. Nature 586, 769-775, doi:10.1038/s41586-020-2786-7 (2020).' 
# zcat /medpop/esp/mesbah/GWAS_CHIP/MPN_summary/MPN_metaGWAS_sumstats.tsv.gz | head -2
# MarkerName	RSID	CHR	POS	REF	ALT	Effect	StdErr	pvalue	MAF
# 1:768448_G_A	rs12562034	1	768448	G	A	-0.0253	0.0615	0.681	0.1076
# Effect allele = ALT allele
# MPN N = 1156774
mpn <- fread("/Volumes/mesbah/gwas/mashr/chr1_22.mpn_hg37_summary.tsv.gz",
             select = c("SNPID", "MPN_Z", "MPN_N"))

## LTL GWAS
# zcat /medpop/esp/mesbah/GWAS_CHIP/ltl_gwas/UKB_telomere_gwas_summarystats.tsv.gz | head -2                                                       
# variant_id	p_value	chromosome	base_pair_location	effec_allele	other_allele	effect_allele_frequency	beta	standard_error
# rs367896724	0.13	1	10177	A	AC	0.599214	0.00452419	0.00296934
ltl <- fread("/Volumes/mesbah/gwas/mashr/chr1_22.ltl_hg37_summary.tsv.gz",
             select = c("SNPID", "LTL_Z", "LTL_N"))

## mCA GWAS
mca <- fread("/Volumes/mesbah/gwas/mashr/chr1_22.expanded_mCA_hg37_summary.tsv.gz",
             select = c("SNPID", "mCA_Z", "mCA_N"))

# save.image(file = "/Volumes/mesbah/gwas/mashr/mpn_ltl_mCA.4mash.rda")
# load("/Volumes/mesbah/gwas/mashr/mpn_ltl_mCA.4mash.rda")

## Overall CHIP GWAS
chip <- fread("/Volumes/mesbah/gwas/mashr/chr1_22.CHIP_hg37_summary.tsv.gz",
             select = c("SNPID", "CHIP_Z", "CHIP_N"))

chip <- subset(chip, chip$SNPID %in% unique(c(mpn$SNPID, ltl$SNPID, mca$SNPID)) )

## DNMT3A CHIP GWAS
dnmt3a <- fread("/Volumes/mesbah/gwas/mashr/chr1_22.DNMT3A_hg37_summary.tsv.gz",
             select = c("SNPID", "DNMT3A_Z", "DNMT3A_N"))

dnmt3a <- subset(dnmt3a, dnmt3a$SNPID %in% unique(c(mpn$SNPID, ltl$SNPID, mca$SNPID)) )


## TET2 CHIP GWAS
tet2 <- fread("/Volumes/mesbah/gwas/mashr/chr1_22.TET2_hg37_summary.tsv.gz",
             select = c("SNPID", "TET2_Z", "TET2_N"))

tet2 <- subset(tet2, tet2$SNPID %in% unique(c(mpn$SNPID, ltl$SNPID, mca$SNPID)) )


## merge
  # TET2 and CHIP
dat <- merge(tet2, chip, by="SNPID",  all.x = T); rm(tet2, chip)

  # DNMT3A
dat <- merge(dat, dnmt3a, by="SNPID",  all.x = T); rm(dnmt3a)

  # LTL
dat <- merge(dat, ltl, by="SNPID",  all.x = T); rm(ltl)


  # mCA
dat <- merge(dat, mca, by="SNPID",  all.x = T); rm(mca)


  # MPN
dat <- merge(dat, mpn, by="SNPID",  all.x = T); rm(mpn)

## Save 
# save(dat, file = "/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/mash/mash_file.chip_dnmt3a_tet2_mpn_ltl_mCA.rda")

# load("/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/mash/mash_file.chip_dnmt3a_tet2_mpn_ltl_mCA.rda")

## correlation
# cor(dat[,c(2,4,6,8,10,12)], use="complete")
# 
# dat$P <- pnorm(abs(dat$CHIP_Z), lower.tail = FALSE)
# 
# dat_p5e7 <- subset(dat, dat$P<5e-7)
# 
# heatmap(cor(dat_p5e7[,c(2,4,6,8,10,12)], use="complete"))
