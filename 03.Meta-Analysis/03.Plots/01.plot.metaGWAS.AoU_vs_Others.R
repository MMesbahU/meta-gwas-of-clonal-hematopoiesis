
######### Read inputs from command line: gwas.gz outpath trait prefix_fig Sample_size minMAF maxMAF
ARGs <- commandArgs(TRUE)
###
library(data.table)

library(qqman)
###
# V1
inputGWAS <- ARGs[1]

outDir <- ARGs[2]

trait <- ARGs[3]

prefix <- ARGs[4]

Sample_size <- as.numeric(ARGs[5])

minMAF <- as.numeric(ARGs[6])

maxMAF <- as.numeric(ARGs[7])


setwd(outDir)

## Any CHIP
meta_gwas <- fread(cmd=paste("gsutil cat ", inputGWAS," | zcat"), header=T)
## varID	EffectAllele	OtherAllele	EAF	BETA	SE	P	q_P	N	Effect_Direction
# names(meta_gwas) <- c("SNP","Chr", "POS", "P")

## Filter by MAF and total N
# N>=700k
cat(Sample_size,"\n")
## 
meta_gwas <- subset(meta_gwas,
		    meta_gwas$N >= Sample_size & 
		    meta_gwas$EAF >= minMAF &
                    meta_gwas$EAF <= maxMAF)

###
library(stringr)

meta_gwas$CHR <- (str_split_fixed(string = meta_gwas$varID, pattern = "[:]", n = 4)[,1])

meta_gwas$CHR <- as.numeric(gsub(pattern="chr", replacement="", x=meta_gwas$CHR))

meta_gwas$POS <- as.numeric(str_split_fixed(string = meta_gwas$varID, pattern = "[:]", n = 4)[,2])

cat("Number of SNVs: ", nrow(meta_gwas),"\n")

head(meta_gwas)


# Autosomes
meta_gwas <- subset(meta_gwas, meta_gwas$CHR %in% c(1:22) ) 

meta_gwas <- meta_gwas[, c(1, 11, 12, 7)]

#####
pdf(paste0("manhattan.metaGWAS.",trait,".",prefix,".pdf"), 
    width = 12, 
    height = 7)

manhattan(x = meta_gwas , chr = "CHR", bp = "POS", p = "P", 
          snp = "varID", suggestiveline=F, 
          col = c("darkblue", "gray60"), 
          main = paste0("GWAS for ",trait), 
          ylim=c(0, round(max(-log10(meta_gwas$P))+5,-1)))

dev.off()

  # QQ
pdf(paste0("qqplot.metaGWAS.",trait,".",prefix,".pdf"), width = 5, height = 6)

gwas_lambda <- round(median(qchisq(1 - meta_gwas$P, 1), na.rm = T)/qchisq(0.5,1), 2)

cat(trait," lambda=", gwas_lambda,"\n")

qq(meta_gwas$P, main=paste0(trait,": lambda = ", gwas_lambda) )

dev.off()

##########
