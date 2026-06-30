
################## Compare CHIP GWAS Loci from 
################## Kar et al. 2022 Nat Gen ########### 
## GCST90102618 (overall CH; https://www.ebi.ac.uk/gwas/studies/GCST90102618)
## GCST90102619 (DNMT3A-CH; https://www.ebi.ac.uk/gwas/studies/GCST90102619) 
## GCST90102620 (TET2-CH; https://www.ebi.ac.uk/gwas/studies/GCST90102620)
library(data.table)
library(liftOver)
library(rtracklayer)
# ch.kr <- fread ("/Volumes/mesbah/gwas/ukb450k/kar.et.a.2022/GCST90102618_buildGRCh37.tsv.gz")
# n snps = 13,392,859
# ch.kr$fdr <- p.adjust(p = ch.kr$p_value, method = "BH")
# ch.kr <- subset(ch.kr, ch.kr$fdr<0.05)
# ch.kr$chromosome <- paste0("chr",ch.kr$chromosome)
# ch.kr.sig5e8 <- subset(ch.kr, ch.kr$p_value<5e-8)
# variant_id	chromosome	base_pair_location	effect_allele	other_allele	effect_allele_frequency	INFO	beta	standard_error	CHISQ_BOLT_LMM	p_value
# rs367896724	1	10177	A	AC	0.600931	0.467935	5.82621e-4	8.32341e-4	0.456309	0.5
ch.kr.sig5e8 <- fread ("gzcat /Volumes/mesbah/gwas/ukb450k/kar.et.a.2022/chip.GCST90102618_buildGRCh37.tsv.gz | awk \'(NR==1){print $0}(NR>1 && $11<5e-8){print $0}\'")
ch.kr.sig5e8$end <- ch.kr.sig5e8$base_pair_location
ch.kr.sig5e8$Var_hg37 <- paste(ch.kr.sig5e8$chromosome, ch.kr.sig5e8$base_pair_location,
                               ch.kr.sig5e8$effect_allele, ch.kr.sig5e8$other_allele, 
                               sep=":")
ch.kr.sig5e8 <- ch.kr.sig5e8[, c(2,3,13,4,5,1,6:12,14)]  
names(ch.kr.sig5e8) <- c("chrom", "start", "end", "ref", "alt", "dbSNP_id", names(ch.kr.sig5e8)[7:14] )
# pchisq(172.6110, 1, lower.tail = F)
  # DNMT3A
dnmt.kr <- fread ("/Volumes/mesbah/gwas/ukb450k/kar.et.a.2022/GCST90102619_buildGRCh37.tsv.gz")
dnmt.kr$fdr <- p.adjust(p = dnmt.kr$p_value, method = "BH")
dnmt.kr <- subset(dnmt.kr, dnmt.kr$fdr<0.05)
dnmt.kr$chromosome <- paste0("chr",dnmt.kr$chromosome)
dnmt.kr.sig5e8 <- subset(dnmt.kr, dnmt.kr$p_value<5e-8)

  # TET2
tet.kr <- fread ("/Volumes/mesbah/gwas/ukb450k/kar.et.a.2022/GCST90102620_buildGRCh37.tsv.gz")
tet.kr$fdr <- p.adjust(p = tet.kr$p_value, method = "BH")
tet.kr <- subset(tet.kr, tet.kr$fdr<0.05)
tet.kr$chromosome <- paste0("chr",tet.kr$chromosome)
tet.kr.sig5e8 <- subset(tet.kr, tet.kr$p_value<5e-8)

# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# library(GenomicRanges)
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# tx_hg19 <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
# unzipped chain file needed
# Load Chain file
# chain = import.chain("~/Documents/Project/Baylor_ARIC_Exomes/hg19ToHg38.over.chain")
  # G-Ranges
# gr.ch.kr.sig5e8 <- as(ch.kr.sig5e8, "GenomicRanges")

# ch.kr.sig5e8.hg38 <- liftOver(gr.ch.kr.sig5e8, chain)

# ch.kr.sig5e8.hg38 <- as.data.frame(ch.kr.sig5e8.hg38)
# ch.kr.sig5e8.hg38$hg38varID <- paste(ch.kr.sig5e8.hg38$seqnames, 
#                                      ch.kr.sig5e8.hg38$start, 
#                                      ch.kr.sig5e8.hg38$ref,
#                                      ch.kr.sig5e8.hg38$alt, 
#                                     sep = ":")

# overlap with GWAS 
# length(unique(ch.kr.sig5e8.hg38$hg38varID))

table(unique(ch.kr.sig5e8.hg38$hg38varID) %in% CHIP$RSID)
# TRUE 
# 932

#### Our GWAS
ch.mmu <- fread("/Volumes/mesbah/gwas/ukb450k/metagwas/results/metaGWAS.CHIP.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.out.gz")
dnmt3a.mmu <- fread("/Volumes/mesbah/gwas/ukb450k/metagwas/results/metaGWAS.DNMT3A.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.out.gz")
tet2.mmu <- fread("/Volumes/mesbah/gwas/ukb450k/metagwas/results/metaGWAS.TET2.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.out.gz")


### Kessler et al. 2022
library(readxl)
rgc.ch <- read_excel("/Users/muddin/Documents/Project/CHIP_GWAS/Manuscript/Draft/Literetures/Regeneron.media-1.xlsx",sheet = "S4_UKB_CHIP_inclusive_index_snp")

rgc.dnmt3a <- read_excel("/Users/muddin/Documents/Project/CHIP_GWAS/Manuscript/Draft/Literetures/Regeneron.media-1.xlsx", sheet = 10)

rgc.tet2 <- read_excel("/Users/muddin/Documents/Project/CHIP_GWAS/Manuscript/Draft/Literetures/Regeneron.media-1.xlsx", sheet = 13)
