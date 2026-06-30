
######### Read inputs from command line
ARGs <- commandArgs(TRUE)

library(data.table)

library(qqman)

# V1
inputGWAS <- ARGs[1]

outQQplot_001 <- ARGs[2]
outQQplot_002 <- ARGs[3]
outQQplot_005 <- ARGs[4]
outQQplot_01 <- ARGs[5]

##
cat(inputGWAS,"\n")

## Any CHIP
# SNPID	EAF	Z	P	N
gwas <- fread(inputGWAS, header=T)
cat("EAF\n")
summary(gwas$EAF)
cat("N\n")
summary(gwas$N)
cat("P\n")
summary(gwas$P)
###
gwas002 <- subset(gwas, gwas$EAF>=0.002 & gwas$EAF<=(1-0.002) )
gwas005 <- subset(gwas, gwas$EAF>=0.005 & gwas$EAF<=(1-0.005) )
gwas01 <- subset(gwas, gwas$EAF>=0.01 & gwas$EAF<=(1-0.01) )
###
# png(outQQplot_001, width = 5, height = 6, units = 'in', res = 300)
# did not work : pdf(outQQplot, width = 5, height = 6)

gwas_lambda.eaf001 <- round(median(qchisq(1 - gwas$P, 1), na.rm = T)/qchisq(0.5,1), 3)

gwas_lambda.eaf002 <- round(median(qchisq(1 - gwas002$P, 1), na.rm = T)/qchisq(0.5,1), 3)

gwas_lambda.eaf005 <- round(median(qchisq(1 - gwas005$P, 1), na.rm = T)/qchisq(0.5,1), 3)

gwas_lambda.eaf01 <- round(median(qchisq(1 - gwas01$P, 1), na.rm = T)/qchisq(0.5,1), 3)

cat("EAF 0.1% lambda=", gwas_lambda.eaf001,"\n")
cat("EAF 0.2% lambda=", gwas_lambda.eaf002,"\n")
cat("EAF 0.5% lambda=", gwas_lambda.eaf005,"\n")
cat("EAF 1% lambda=", gwas_lambda.eaf01,"\n")

png(outQQplot_001, width = 5, height = 6, units = 'in', res = 300)
qq(gwas$P, main=paste0("EAF 0.1% lambda = ", round(gwas_lambda.eaf001, 2) ) )
dev.off()

png(outQQplot_002, width = 5, height = 6, units = 'in', res = 300)
qq(gwas002$P, main=paste0("EAF 0.2% lambda = ", round(gwas_lambda.eaf002, 2) ) )
dev.off()

png(outQQplot_005, width = 5, height = 6, units = 'in', res = 300)
qq(gwas005$P, main=paste0("EAF 0.5% lambda = ", round(gwas_lambda.eaf005, 2) ) )
dev.off()

png(outQQplot_01, width = 5, height = 6, units = 'in', res = 300)
qq(gwas01$P, main=paste0("EAF 1% lambda = ", round(gwas_lambda.eaf01, 2)) )
dev.off()

