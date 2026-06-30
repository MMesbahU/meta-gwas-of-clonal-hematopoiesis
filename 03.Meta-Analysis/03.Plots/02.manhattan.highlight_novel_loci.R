###### Manhattan plot
library(data.table)
# library(qqman)
library(calibrate)
# source("~/Documents/Project/CHIP_GWAS/highlight_new_oldLoci.manhattan.R")
source("/Volumes/medpop_esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/03.Plots/highlight_new_oldLoci.manhattan.R")
# 
# chip <- fread("~/Documents/Project/CHIP_GWAS/Manhattan_plots/sumstats/chr1_22.lifted_hg37.GWAMA.meta_CHIP.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.tsv.gz", 
#               skip=1, header=F)
CHIP <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/multi_ancestry_summary/lifted_hg37.metaGWAS.CHIP.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.hg37_dbSNP.eaf001_min2Studies.tsv.gz", 
              select = c("MarkerID","CHR",	"POS"	,"P"))
names(CHIP) <- c("SNP","CHR", "POS", "P")
nrow(CHIP)
# n SNPs = 23802037

DNMT3A <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/multi_ancestry_summary/lifted_hg37.metaGWAS.DNMT3A.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.hg37_dbSNP.eaf001_min2Studies.tsv.gz", 
              select = c("MarkerID","CHR",	"POS"	,"P"))
names(DNMT3A) <- c("SNP","CHR", "POS", "P")
nrow(DNMT3A)
#n snps = 23883996
TET2 <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/multi_ancestry_summary/lifted_hg37.metaGWAS.TET2.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.hg37_dbSNP.eaf001_min2Studies.tsv.gz", 
              select = c("MarkerID","CHR",	"POS"	,"P"))
names(TET2) <- c("SNP","CHR", "POS", "P")
nrow(TET2)
# n SNPs = 23902838
## DNMT3A
# chip <- fread("~/Documents/Project/CHIP_GWAS/Manhattan_plots/sumstats/chr1_22.lifted_hg37.GWAMA.meta_DNMT3A.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.tsv.gz", 
#               skip=1, header=F)
## TET2
# chip <- fread("~/Documents/Project/CHIP_GWAS/Manhattan_plots/sumstats/chr1_22.lifted_hg37.GWAMA.meta_TET2.allsamples_topmed_ukbb_mgbb_BioVU.eaf001_min2Studies.tsv.gz", 
#               skip=1, header=F)
# CHIP <- chip [,c(5, 1, 2, 12)]; rm(chip)

# Locus overlap 
overlap_chip <- fread("~/Documents/Project/CHIP_GWAS/rerun/overlap/summary_p5e8.CHIP.csv")
table(overlap_chip$Gene, overlap_chip$Final)
# overlap_chip <- merge(overlap_chip, topHits, by.x="MarkerID", by.y="SNP")
# CHIP$SNP[(CHIP$SNP)=="chr12:26422298:GAAT:G"] <- "rs143848555"
trait="Overall CHIP"
CHIPSigP_SNP <- CHIP$SNP[CHIP$P<5e-8]

CHIPSigP_SNP_new <- CHIPSigP_SNP [CHIPSigP_SNP%in%overlap_chip$MarkerID[overlap_chip$Final=="Y"]]  

CHIPSigP_SNP_old <- CHIPSigP_SNP [CHIPSigP_SNP%in%overlap_chip$MarkerID[overlap_chip$Final=="N"]]  

## DNMT3A
# overlap_dnmt3a <- fread("~/Documents/Project/CHIP_GWAS/rerun/overlap/summary_p5e8.DNMT3A.csv")
# table(overlap_dnmt3a$Gene, overlap_dnmt3a$Final)
# trait="DNMT3A"
# SigP_SNP <- DNMT3A$SNP[DNMT3A$P<=5e-8]
# SigP_SNP_new <- SigP_SNP [SigP_SNP %in% overlap_dnmt3a$MarkerID[overlap_dnmt3a$Final=="Y"]]  
# SigP_SNP_old <- SigP_SNP [SigP_SNP %in% overlap_dnmt3a$MarkerID[overlap_dnmt3a$Final=="N"]]  

# TET2
trait="TET2"
overlap_tet2 <- fread("~/Documents/Project/CHIP_GWAS/rerun/overlap/summary_p5e8.TET2.csv")
table(overlap_tet2$Gene, overlap_tet2$Final)
trait="TET2"
SigP_SNP <- TET2$SNP[TET2$P<=5e-8]
SigP_SNP_new <- SigP_SNP [SigP_SNP %in% overlap_tet2$MarkerID[overlap_tet2$Final=="Y"]]  
SigP_SNP_old <- SigP_SNP [SigP_SNP %in% overlap_tet2$MarkerID[overlap_tet2$Final=="N"]]  

d <- TET2

d <- d[order(d$P), ]

d <- subset(d, !duplicated(d$SNP))

## 
# topHits = subset(d, P <= 5e-8) 
# topHits <- topHits[order(topHits$P), ]
# topSNPs <- NULL
# for (i in unique(topHits$CHR)) {
#   chrSNPs <- topHits[topHits$CHR == i, ]
#   topSNPs <- rbind(topSNPs, chrSNPs[1, ])
# }

# topHits = subset(overlap_chip, P <= 5e-8) 
# topHits <- topHits[order(topHits$P), ]
# topSNPs <- NULL
# for (i in unique(topHits$Locus)) {
#   chrSNPs <- topHits[topHits$Locus == i, ]
#   topSNPs <- rbind(topSNPs, chrSNPs[1, ])
# }
topHits = subset(overlap_dnmt3a, P <= 5e-8) 
topHits <- topHits[order(topHits$P), ]
topSNPs <- NULL
for (i in unique(topHits$Locus)) {
  chrSNPs <- topHits[topHits$Locus == i, ]
  topSNPs <- rbind(topSNPs, chrSNPs[1, ])
}
# d$SNP[(d$SNP)=="chr12:26422298:GAAT:G"] <- "rs143848555"

# png(paste0("~/Documents/Project/CHIP_GWAS/Manhattan_plots/manhattan.",trait,".gwas.meta.maf001_n2_SigP167e8.highlight.png"), 
#     width = 12, height = 7,units = 'in', res = 300 )
png(paste0("~/Documents/Project/CHIP_GWAS/rerun/Figures/manhattan.",trait,".gwas.meta.maf001_n2_SigP5e8.highlight.png"), 
    width = 12, height = 7,units = 'in', res = 300 )
# annotateTop=T, annotatePval = 5e-8,
highlight_manhattan(x = d , chr = "CHR", bp = "POS", p = "P", 
                    snp = "SNP", suggestiveline=F, 
                    genomewideline = -log10(5e-08),
                    highlight_new = SigP_SNP_new, 
                    col_new_highlight = "red",
                    highlight_old = SigP_SNP_old, 
                    col_old_highlight = "blue", 
                    col = c("darkblue", "gray60"), 
                    main = paste0(trait), 
                    ylim=c(0, round(max(-log10(d$P))+5,-1)))

dev.off()
  ## QQ-plot
# png(paste0("~/Documents/Project/CHIP_GWAS/Manhattan_plots/qqplot.",trait,".gwas.meta.maf001_n2_SigP167e8.png"), width = 5, height = 6, units = 'in', res = 300)
png(paste0("~/Documents/Project/CHIP_GWAS/rerun/Figures/qqplot.",trait,".gwas.meta.maf001_n2_SigP5e8.png"), 
    width = 5, height = 6, units = 'in', res = 300)
gwas_lambda <- round(median(qchisq(1 - d$P, 1), na.rm = T)/qchisq(0.5,1), 2)
# chipgwas_lambda <- round(median(qchisq(1 - CHIP$P, 1), na.rm = T)/qchisq(0.5,1), 2)
# dnmtgwas_lambda <- round(median(qchisq(1 - DNMT3A$P, 1), na.rm = T)/qchisq(0.5,1), 2)
# tet2gwas_lambda <- round(median(qchisq(1 - TET2$P, 1), na.rm = T)/qchisq(0.5,1), 2)

cat(trait," lambda=", gwas_lambda,"\n")
# qqman::qq(d$P, main=paste0(trait,": lambda = ", gwas_lambda) )
# qqman::qq(d$P, main=expression(~ lambda ~ "= 1.05") ) # CHIP
# qqman::qq(d$P, main=expression(~ lambda ~ "= 1.04") ) # DNMT3A
qqman::qq(d$P, main=expression(~ lambda ~ "= 1.05") ) # TET2
dev.off()
