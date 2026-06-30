
######### Read inputs from command line
ARGs <- commandArgs(TRUE)

library(data.table)
library(qqman)

# V1
inputGWAS <- ARGs[1]

outDir <- ARGs[2]

trait <- ARGs[3]

prefix <- ARGs[4]

setwd(outDir)

## Any CHIP
gwas <- fread(inputGWAS, header=T)

names(gwas) <- c("SNP","Chr", "POS", "P")

## Convert to numeric 
gwas$CHR <- gsub(pattern="chr", replacement="", x=gwas$Chr)
# gwas$CHR[gwas$CHR == "X"] <- 23

# Only Autosomes
gwas <- subset(gwas, gwas$CHR %in% c(1:22) ) 

gwas$CHR <- as.numeric(gwas$CHR)

# gwas$POS <- as.numeric(gwas$POS)

png(paste0("manhattan.",trait,".",prefix,".gwas.meta.png"), width = 12, height = 7,units = 'in', res = 300 )
manhattan(x = gwas , chr = "CHR", bp = "POS", p = "P", 
          snp = "SNP", suggestiveline=F, 
          col = c("darkblue", "gray60"), 
          main = paste0("GWAS for ",trait), 
          ylim=c(0, round(max(-log10(gwas$P))+5,-1)))
dev.off()

  # QQ
png(paste0("qqplot.",trait,".",prefix,".gwas.meta.png"), width = 5, height = 6, units = 'in', res = 300)
gwas_lambda <- round(median(qchisq(1 - gwas$P, 1), na.rm = T)/qchisq(0.5,1), 2)
cat(trait," lambda=", gwas_lambda,"\n")
qq(gwas$P, main=paste0(trait,": lambda = ", gwas_lambda) )
dev.off()


