##
# UpSet plot of CHIP, MPN, mCA, LTL GWAS

library(UpSetR)
library(data.table)

# chip <- read.csv("~/Documents/Project/CHIP_GWAS/Manuscript/Draft/CHIP2.0_Manuscript/Overlaps/p5e8.chr1_22.lifted_hg37.GWAMA.meta_CHIP.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.csv", stringsAsFactors = F)
# chip.p167e8 <- subset(chip, chip$Pvalue<=1.67e-8)
# chip.p167e8$chr_pos <- paste(chip.p167e8$CHR, chip.p167e8$POS, sep = ":")
# 
# dnmt3a <- read.csv("~/Documents/Project/CHIP_GWAS/Manuscript/Draft/CHIP2.0_Manuscript/Overlaps/p5e8.chr1_22.lifted_hg37.GWAMA.meta_DNMT3A.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.csv", stringsAsFactors = F)
# dnmt3a.p167e8 <- subset(dnmt3a, dnmt3a$Pvalue<=1.67e-8)
# dnmt3a.p167e8$chr_pos <- paste(dnmt3a.p167e8$CHR, dnmt3a.p167e8$POS, sep = ":")
# 
# tet2 <- read.csv("~/Documents/Project/CHIP_GWAS/Manuscript/Draft/CHIP2.0_Manuscript/Overlaps/p5e8.chr1_22.lifted_hg37.GWAMA.meta_TET2.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.csv", stringsAsFactors = F)
# tet2.p167e8 <- subset(tet2, tet2$Pvalue<=1.67e-8)
# tet2.p167e8$chr_pos <- paste(tet2.p167e8$CHR, tet2.p167e8$POS, sep = ":")

mu.ch <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/summary_p5e8.CHIP.tsv")
mu.ch$chr_pos <- paste(mu.ch$CHR,mu.ch$POS, sep = ":")

mu.dnmt3a <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/summary_p5e8.DNMT3A.tsv")
mu.dnmt3a$chr_pos <- paste(mu.dnmt3a$CHR,mu.dnmt3a$POS, sep = ":")

mu.tet2 <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/summary_p5e8.TET2.tsv")
mu.tet2$chr_pos <- paste(mu.tet2$CHR,mu.tet2$POS, sep = ":")

mpn <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/p5e8.MPN_metaGWAS_sumstats.tsv", stringsAsFactors = F)
mpn$chr_pos <- paste(mpn$CHR, mpn$POS, sep = ":")

# ltl <- read.csv("~/Documents/Project/CHIP_GWAS/Manuscript/Draft/CHIP2.0_Manuscript/Overlaps/p5e8.UKB_telomere_gwas_summarystats.csv", stringsAsFactors = F)
ltl <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/p5e8.UKB_telomere_gwas_summarystats.tsv", stringsAsFactors = F)
ltl$chr_pos <- paste(ltl$chromosome, ltl$base_pair_location, sep = ":")

# mca <- read.csv("~/Documents/Project/CHIP_GWAS/Manuscript/Draft/CHIP2.0_Manuscript/Overlaps/p5e8.Maryam.expanded_mCA_summary_stats.N444199.csv", stringsAsFactors = F)
mca <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/p5e8.Maryam.expanded_mCA_summary_stats.N444199.tsv")
mca$chr_pos <- paste(mca$`#chrom`, mca$pos, sep = ":")

# ## merge 
# chip_ltl <- merge(chip, ltl, by.x="RSID", by.y = "variant_id")
# 
# 
# table(unique(chip.p167e8$CHR[chip.p167e8$RSID %in% ltl$variant_id]))
# table(unique(dnmt3a.p167e8$CHR[dnmt3a.p167e8$RSID %in% ltl$variant_id]))
# table(unique(tet2.p167e8$CHR[tet2.p167e8$RSID %in% ltl$variant_id]))
# 
# table(unique(chip.p167e8$CHR[chip.p167e8$RSID %in% mpn$RSID]))
# table(unique(dnmt3a.p167e8$CHR[dnmt3a.p167e8$RSID %in% mpn$RSID]))
# table(unique(tet2.p167e8$CHR[tet2.p167e8$RSID %in% mpn$RSID]))
# 
# table(unique(chip.p167e8$CHR[chip.p167e8$RSID %in% mca$rsid]))
# table(unique(dnmt3a.p167e8$CHR[dnmt3a.p167e8$RSID %in% mca$rsid]))
# table(unique(tet2.p167e8$CHR[tet2.p167e8$RSID %in% mca$rsid]))
# 
# table(chip.p167e8$CHR[chip.p167e8$RSID %in% dnmt3a.p167e8$RSID])
# table(chip.p167e8$CHR[chip.p167e8$RSID %in% tet2.p167e8$RSID])
# table(chip.p167e8$CHR[chip.p167e8$RSID %in% mca$rsid])
# table(chip.p167e8$CHR[chip.p167e8$RSID %in% mpn$RSID])
# 
# table(dnmt3a$CHR[dnmt3a$RSID %in% ltl$variant_id])
# table(dnmt3a$CHR[dnmt3a$RSID %in% mpn$RSID])

## chromosomes
# chroms <- c(unique(c(mu.ch$CHR,mu.dnmt3a$CHR,mu.tet2$CHR)),1,6)
chroms <- c(unique(c(mu.ch$CHR,mu.dnmt3a$CHR,mu.tet2$CHR)))

chip_overlap <- as.data.frame(matrix(NA, nrow = length(chroms),ncol = 7))
names(chip_overlap) <- c("Locus", "Overall CHIP", "DNMT3A CHIP", "TET2 CHIP", "MPN", "LTL", "Expanded mCA")

chip_overlap$Locus <- chroms 
chip_overlap$`Overall CHIP` <- c(ifelse(chip_overlap$Locus[1:13] %in% unique(chip.p167e8$CHR), 1, 0),0,0)

chip_overlap$`DNMT3A CHIP` <- c(ifelse(chip_overlap$Locus[1:13] %in% unique(dnmt3a.p167e8$CHR), 1, 0),0,0)

chip_overlap$`TET2 CHIP` <- c(0,0,0,1,0,0,0,0,0,0,1,1,1,1,1)

chip_overlap$LTL <- c(1,1,0,1,1,1,0,0,1,0,1,0,0,0,0)

chip_overlap$MPN <- c(0,0,1,1,0,1,0,0,0,0,0,0,0,0,0)

chip_overlap$`Expanded mCA` <- c(0,0,0,1,1,1,1,0,1,0,1,0,0,0,0)

### chr_pos overlap
# all_chr_pos <- unique(c(chip.p167e8$chr_pos, dnmt3a.p167e8$chr_pos, tet2.p167e8$chr_pos, mpn$chr_pos, ltl$chr_pos, mca$chr_pos))

all_chip_chr_pos <- unique(c(mu.ch$chr_pos, mu.dnmt3a$chr_pos, mu.tet2$chr_pos))

all_chr_pos_overlap <- as.data.frame(matrix(NA, nrow = length(all_chip_chr_pos),ncol = 7))
names(all_chr_pos_overlap) <- c("CHR_POS", "Overall CHIP", "DNMT3A CHIP", "TET2 CHIP", "MPN", "LTL", "Expanded mCA")

all_chr_pos_overlap$CHR_POS <- all_chip_chr_pos

all_chr_pos_overlap$`Overall CHIP` <- ifelse(all_chr_pos_overlap$CHR_POS %in% mu.ch$chr_pos, 1,0)

all_chr_pos_overlap$`DNMT3A CHIP` <- ifelse(all_chr_pos_overlap$CHR_POS %in% mu.ch$chr_pos, 1,0)

all_chr_pos_overlap$`TET2 CHIP` <- ifelse(all_chr_pos_overlap$CHR_POS %in% mu.tet2$chr_pos, 1,0)

all_chr_pos_overlap$MPN <- ifelse(all_chr_pos_overlap$CHR_POS %in% mpn$chr_pos, 1,0)

all_chr_pos_overlap$LTL <- ifelse(all_chr_pos_overlap$CHR_POS %in% ltl$chr_pos, 1,0)

all_chr_pos_overlap$`Expanded mCA` <- ifelse(all_chr_pos_overlap$CHR_POS %in% mca$chr_pos, 1,0)

## Format data for upset plot

# png("~/Documents/Project/CHIP_GWAS/Manuscript/Draft/CHIP2.0_Manuscript/Overlaps/upset_fig2d.png", 
#     width=7, height=5, units= "in", res=300, pointsize = 4)
# par(mfrow=c(1,1), mar= c(5, 4, 4, 2))
# upset(chip_overlap, sets = c("Overall CHIP", "DNMT3A CHIP", "TET2 CHIP", "MPN", "LTL", "Expanded mCA"),
#       order.by ="freq", decreasing = F,
#       mb.ratio = c(0.7, 0.3),
#       sets.bar.color=c("red", "cornflowerblue", "seagreen", "gray45", "black", "orange"))
# dev.off()
# pdf("~/Documents/Project/CHIP_GWAS/Manuscript/Draft/CHIP2.0_Manuscript/Overlaps/upset_fig2d.pdf",
#     width=7, height=5, bg="transparent")
pdf("~/Documents/Project/CHIP_GWAS/rerun/overlap/upset_fig2d.pdf",
    width=7, height=5, bg="transparent")

upset(all_chr_pos_overlap, sets = c("Overall CHIP", 
                                    "DNMT3A CHIP", 
                                    "TET2 CHIP", 
                                    "MPN", 
                                    "LTL", 
                                    "Expanded mCA"),
      order.by ="freq", decreasing = T, 
      mb.ratio = c(0.7, 0.3))

dev.off()

## Save overlap data
# save.image(file = "~/Documents/Project/CHIP_GWAS/Manuscript/Draft/CHIP2.0_Manuscript/Overlaps/overlap_chip_ltl_mpn_mCA.rda")

################################ CHIP/DNMT3A/TET2 Loci Overlap
### Loci Overlap
mu.ch <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/summary_p5e8.CHIP.tsv")
mu.ch$chr_pos <- paste(mu.ch$CHR,mu.ch$POS, sep = ":")

mu.dnmt3a <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/summary_p5e8.DNMT3A.tsv")
mu.dnmt3a$chr_pos <- paste(mu.dnmt3a$CHR,mu.dnmt3a$POS, sep = ":")

mu.tet2 <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/summary_p5e8.TET2.tsv")
mu.tet2$chr_pos <- paste(mu.tet2$CHR,mu.tet2$POS, sep = ":")

########## Kar et al. 2022
kar.ch <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/kar.chip.GCST90102618_buildGRCh37.p5e8.tsv")
kar.ch$chr_pos <- paste(kar.ch$chromosome,kar.ch$base_pair_location, sep = ":")

kar.dnmt3a <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/kar.dnmt3a.GCST90102619_buildGRCh37.p5e8.tsv")
kar.dnmt3a$chr_pos <- paste(kar.dnmt3a$chromosome,kar.dnmt3a$base_pair_location, sep = ":")

kar.tet2 <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/overlap/kar.tet2.GCST90102620_buildGRCh37.p5e8.tsv")
kar.tet2$chr_pos <- paste(kar.tet2$chromosome,kar.tet2$base_pair_location, sep = ":")

#################
### Kessler et al. 2022
library(readxl)
rgc.ch <- read_excel("/Users/muddin/Documents/Project/CHIP_GWAS/Manuscript/Draft/Literetures/Regeneron.media-1.xlsx",sheet = "S4_UKB_CHIP_inclusive_index_snp")

rgc.dnmt3a <- read_excel("/Users/muddin/Documents/Project/CHIP_GWAS/Manuscript/Draft/Literetures/Regeneron.media-1.xlsx", sheet = 10)

rgc.tet2 <- read_excel("/Users/muddin/Documents/Project/CHIP_GWAS/Manuscript/Draft/Literetures/Regeneron.media-1.xlsx", sheet = 13)
#################