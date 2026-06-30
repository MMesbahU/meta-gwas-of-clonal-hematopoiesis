### Create genomic intervals
# https://stackoverflow.com/questions/39458808/break-region-into-smaller-regions-based-on-cutoff
library(GenomicRanges)
hg38 <- fread("/Volumes/medpop_esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/Homo_sapiens_assembly38.fasta.fai", nrows = 22)
hg38$Start <- 1
hg38$End <- hg38$V2
hg38.1_10 <-  hg38[1:10,]
hg38.11_22 <- hg38[11:22,]

rngs.hg38.1_10 = GRanges(hg38.1_10$V1, IRanges(hg38.1_10$Start, hg38.1_10$End))
tiles.hg38.1_10 = tile(rngs.hg38.1_10, n=10)
hg38_range.1_10 <- as(unlist(tiles.hg38.1_10), "data.frame")

rngs.hg38.11_22 = GRanges(hg38.11_22$V1, IRanges(hg38.11_22$Start, hg38.11_22$End))
tiles.hg38.11_22 = tile(rngs.hg38.11_22, n=5)
hg38_range.11_22 <- as(unlist(tiles.hg38.11_22), "data.frame")

hg38_range <- as.data.frame(rbind(hg38_range.1_10, hg38_range.11_22))
hg38_range$chrom_num <- as.numeric(gsub(pattern = "chr", replacement = "",x = hg38_range$seqnames))
hg38_range$chrom_interval <- paste(hg38_range$start, hg38_range$end, 
                                   sep = "-")
hg38_range <- hg38_range[,c(6,7,1:4)]

fwrite(hg38_range, "/Volumes/medpop_esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg38_chrom_interval.tsv", 
				          row.names = F, col.names = F, quote = F, sep="\t")

###### hg37/hg19:
hg19 <- fread("/Volumes/medpop_esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/Homo_sapiens_assembly19.fasta.fai", nrows = 22)
hg19$Start <- 1
hg19$End <- hg19$V2
hg19.1_10 <-  hg19[1:10,]
hg19.11_22 <- hg19[11:22,]

rngs.hg19.1_10 = GRanges(hg19.1_10$V1, IRanges(hg19.1_10$Start, hg19.1_10$End))
tiles.hg19.1_10 = tile(rngs.hg19.1_10, n=10)
hg19_range.1_10 <- as(unlist(tiles.hg19.1_10), "data.frame")

rngs.hg19.11_22 = GRanges(hg19.11_22$V1, IRanges(hg19.11_22$Start, hg19.11_22$End))
tiles.hg19.11_22 = tile(rngs.hg19.11_22, n=5)
hg19_range.11_22 <- as(unlist(tiles.hg19.11_22), "data.frame")

hg19_range <- as.data.frame(rbind(hg19_range.1_10, hg19_range.11_22))
hg19_range$chrom_num <- as.numeric(gsub(pattern = "chr", replacement = "",x = hg19_range$seqnames))
hg19_range$chrom_interval <- paste(hg19_range$start, hg19_range$end, 
                                   sep = "-")
hg19_range <- hg19_range[,c(6,7,1:4)]

fwrite(hg19_range, "/Volumes/medpop_esp2/mesbah/tools/CHIP_metaAnalysis/02.GWAS/1.regenie_GWAS/hg19_chrom_interval.tsv", 
       row.names = F, col.names = F, quote = F, sep="\t")

####################################

