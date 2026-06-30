
######### Read inputs from command line
ARGs <- commandArgs(TRUE)

library(data.table)
library(qqman)

# V1
inputGWAS <- ARGs[1]

outDir <- ARGs[2]

trait <- ARGs[3]

setwd(outDir)

## Any CHIP
gwas <- fread(inputGWAS, header=T, select = c("RSID","CHR","POS", "Pvalue"))

names(gwas) <- c("SNP","Chr", "POS", "P")

## Convert to numeric 
gwas$CHR <- gsub(pattern="chr", replacement="", x=gwas$Chr)
# gwas$CHR[gwas$CHR == "X"] <- 23

# Only Autosomes
# gwas <- subset(gwas, gwas$CHR %in% c(1:22) ) 

gwas$CHR <- as.numeric(gwas$CHR)

# gwas$POS <- as.numeric(gwas$POS)
	## Highlght Significant P
SigP <- gwas$SNP[gwas$P<=1.67e-8]

png(paste0("manhattan.",trait,".gwas.meta.maf001_n2_SigP167e8.highlight.png"), width = 12, height = 7,units = 'in', res = 300 )
manhattan(x = gwas , chr = "CHR", bp = "POS", p = "P", 
          snp = "SNP", suggestiveline=F, 
	  genomewideline = -log10(1.67e-08),
	  highlight = SigP,         
	  col = c("darkblue", "gray60"), 
          main = paste0("GWAS for ",trait), 
          ylim=c(0, round(max(-log10(gwas$P))+5,-1)))
dev.off()

  # QQ
png(paste0("qqplot.",trait,".gwas.meta.maf001.png"), width = 5, height = 6, units = 'in', res = 300)
gwas_lambda <- round(median(qchisq(1 - gwas$P, 1), na.rm = T)/qchisq(0.5,1), 2)
cat(trait," lambda=", gwas_lambda,"\n")
qq(gwas$P, main=paste0(trait,": lambda = ", gwas_lambda) )
dev.off()
