######### Read inputs from command line
ARGs <- commandArgs(TRUE)

library(data.table)

inputGWAS <- ARGs[1]

# outFDR5 <- ARGs[2]

## Any CHIP
gwas <- fread(inputGWAS, select=c('SNP', 'p'))

cat("GWAS File:", inputGWAS,"\n") 
cat("n SNPs in the File:", nrow(gwas),"\n")
gwas$FDR <- p.adjust(p=gwas$p, method="BH")

cat("GWAS Sig (N SNPs with P<=5e-8) =")
print(table(gwas$p<=5e-8))

cat("FDR sig SNPs (N SNPs with FDR<=0.05)=")
print(table(gwas$FDR<=0.05))

cat("max P with FDR<=0.05=")
print(max(gwas[gwas$FDR<=0.05, 2]))

