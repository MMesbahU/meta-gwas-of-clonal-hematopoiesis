######### select SNPs for null model
## GSA array SNPs to extract from TOPMed
gsa_snps <- fread("GSA-24v3-0_A2.csv.gz", 
                  skip=7, header = T, sep=",")

noduplicate_gsa_snps <- subset(gsa_snps[, c(10, 11)], 
                               gsa_snps$MapInfo>0 & 
                                 gsa_snps$Chr %in% c(seq(1:22)))

names(noduplicate_gsa_snps) <- c("CHROM", "POS")

noduplicate_gsa_snps$CHROM <- paste0("chr", 
                                     noduplicate_gsa_snps$CHROM)

noduplicate_gsa_snps$SNP <- paste(noduplicate_gsa_snps$CHROM, 
                                  noduplicate_gsa_snps$POS, 
                                  sep=":")


noduplicate_gsa_snps <- subset(noduplicate_gsa_snps, 
                               !duplicated(noduplicate_gsa_snps$SNP))

# fwrite(noduplicate_gsa_snps, 
#        "chr1_22.hg38.GSA_24v3_0_A2.tsv", 
#        row.names = F, col.names = T, sep="\t")
