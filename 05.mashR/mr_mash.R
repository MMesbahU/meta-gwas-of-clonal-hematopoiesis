
###### MASH analysis 
######
library(mashr)
library(data.table)
library(stringr)
# znew=read.table("~/Dropbox/summaryZ.mpn_chip_dnmt3a_tet2.tsv.gz",sep="\t",header=T)
# /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_N900k/mash/ch_mash_noNA.rda

load("/Volumes/mesbah-1/projects/Meta_GWAS/MetaGWAS_N900k/mash/ch_mash_noNA.rda")

znew <- ch_mash_noNA

# znew$chr <- str_split_fixed(string = znew$SNPID, pattern = ":", n = 4)[, 1]

# znew$bp <- as.numeric(str_split_fixed(string = znew$SNPID, pattern = ":", n = 4)[,2])

a <- do.call(rbind, strsplit(znew$SNPID, ':'))

chr <- a[,1]

# a2=do.call(rbind, strsplit(chr, 'chr'))
# chr=a2[,2]

bp <- do.call(rbind, strsplit(a[,2], '\\|'))[,1]

znew$a.chr <- chr

znew$a.pos <- as.numeric(as.character(bp))

## LD chunk
bed <- read.table("/Volumes/mesbah-1/projects/Meta_GWAS/MetaGWAS_N900k/mash/ld_chunk.bed")

head(bed)

bed$V1 <- do.call(rbind, strsplit(bed$V1, 'chr'))[,2]

## Get Maximum of Z for each trait
maxes <- apply(znew[,c(2,4,6,8,10,12)],1,function(x){max(abs(x))})

znew$maxes <- maxes

# Get maximum for each ld block
max_block <- data.frame(matrix(ncol = ncol(znew), nrow = nrow(bed)))

colnames(max_block) <- (colnames(znew))

for(i in 1:nrow(bed)){
  chr <- bed[i,1]
  start <- bed[i,2]
  stop <- bed[i,3]
  in_chrom <- znew[znew$a.chr==chr,]##extract those of the same
  goodguys <- in_chrom[in_chrom$a.pos > start & in_chrom$a.pos < stop ,]
  if(nrow(goodguys)>0) {
    z.max <- which.max(goodguys[,"maxes"])
    z_good <- goodguys[z.max,]
    z_good$a.chr <- as.character(z_good$a.chr)
  } else {
    z_good <- rep(0,ncol(max_block))
  }
  #z_good=data.table(z_good,stringsAsFactors = F)
  
  #z_good$=rownames(goodguys)[z.max]
  max_block[i,] <- z_good
  
  print(i)
  
}

max_block <- na.omit(max_block)

# saveRDS(max_block,"/Volumes/mesbah-1/projects/Meta_GWAS/MetaGWAS_N900k/mash/chip_mash_max_block.oct2024.rds")

max_block <- readRDS("/Volumes/mesbah-1/projects/Meta_GWAS/MetaGWAS_N900k/mash/chip_mash_max_block.oct2024.rds")

###############
# source('~/Dropbox/mash_materials/flashscript.R')
zmash <- as.matrix(znew[c(2,4,6,8,10,12)])

rownames(zmash) <- znew$SNPID

#zmash=readRDS("~/Dropbox/cadproject/zmash_mat.rds")
# identify a random subset of 40000 tests
set.seed(20241031)

random.subset <- sample(1:nrow(zmash),40000)

data.temp <- mash_set_data(zmash[random.subset,],alpha = 1)

##correlation matrix of the errors
Vhat <- estimate_null_correlation_simple(data.temp) 

# saveRDS(Vhat,"/Volumes/mesbah-1/projects/Meta_GWAS/MetaGWAS_N900k/mashVhat.oct2024.rds")

##
library("lattice")

clrs = colorRampPalette((c("#D73027","#FC8D59","#FEE090","#FFFFBF", "#E0F3F8","#91BFDB","#4575B4")))(64)

print(levelplot(Vhat,
                col.regions = clrs,
                xlab = "",
                ylab = "",
                colorkey = TRUE,
                main="VHAT"))

rm(data.temp)

##use for E(Z)model
data.random <- mash_set_data(zmash[random.subset,],alpha = 1,V=Vhat)

max_block <- readRDS("/Volumes/mesbah-1/projects/Meta_GWAS/MetaGWAS_N900k/mash/chip_mash_max_block.oct2024.rds")

zmax <- apply(max_block[,c(2,4,6,8,10,12)],
              2,function(x){as.numeric(x)})

rownames(zmax) <- max_block$SNPID

data.strong <- mash_set_data(zmax,
                             alpha = 1,
                             V=Vhat)

U.pca <- cov_pca(data.strong, 3)

U.flash <- cov_flash(data.strong,
                     remove_singleton = T)
#
gc()
#, non_canonical = TRUE)
X.center <- apply(data.strong$Bhat, 2, function(x) x - mean(x))

U.ed <- cov_ed(data.strong, 
               c(U.pca,U.flash,
                 list("XX" = t(X.center) %*% X.center / nrow(X.center))))

# saveRDS(U.ed,"/Volumes/mesbah-1/projects/Meta_GWAS/MetaGWAS_N900k/mash/EDcov.oct2024.Rds")
### Run canonical model
U.c <- cov_canonical(data.random)

m <- mash(data.random, Ulist = c(U.ed, U.c), outputlevel = 1)

# saveRDS(m, "/Volumes/mesbah-1/projects/Meta_GWAS/MetaGWAS_N900k/mash/mfitCHIP.oct2024.rds")

# saveRDS(zmash, "/Volumes/mesbah-1/projects/Meta_GWAS/MetaGWAS_N900k/mash/zmash.oct2024.rds")

# m <- readRDS("/Volumes/mesbah-1/projects/Meta_GWAS/MetaGWAS_N900k/mash/mfitCHIP.oct2024.rds")

# zmash <- readRDS("/Volumes/mesbah-1/projects/Meta_GWAS/MetaGWAS_N900k/mash/zmash.oct2024.rds")

k <- length(m$fitted_g$Ulist)

l <- length(m$fitted_g$grid)

pimat <- matrix(m$fitted_g$pi[-1],
                nrow=l,
                byrow=T)

colnames(pimat) <- names(m$fitted_g$pi)[2:(k+1)]

barplot(colSums(pimat),las=2)

############
library("lattice")

for(i in 1:7){
  z.num <- as.matrix(cov2cor(m$fitted_g$Ulist[[i]]))
  colnames(z.num) <- colnames(zmash)[c(1:6)]
  rownames(z.num) <- NULL
  clrs = colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#FFFFBF", "#E0F3F8","#91BFDB","#4575B4")))(64)
  z.num[lower.tri(z.num)] <- NA
  print(levelplot(z.num,col.regions = clrs,xlab = "",ylab = "",colorkey = TRUE,main=paste0(names(m$fitted_g$Ulist)[[i]])))
}

head(z.num)

#########
# Posterior
for(i in 0:2){
  
  start <- i*2e6+1
  
  stop <- (i+1)*2e6
  
  print(c(start, stop))
  
  # library("mashr")
  
  mash.data <- mash_set_data(zmash[start:stop,],
                             V = Vhat,
                             alpha = 1)
  
  p <- mash_compute_posterior_matrices(m$fitted_g, 
                                       mash.data, 
                                       algorithm.version = "Rcpp")
  
  saveRDS(p,
          file = paste0("/Volumes/mesbah-1/projects/Meta_GWAS/MetaGWAS_N900k/mash/chip_mash_batch",i,".oct2024.rds"))
}

## 
i <- 2

start <- i*2e6+1

stop <- nrow(zmash)

print(c(start, stop))

# library("mashr")
mash.data <- mash_set_data(zmash[start:stop,],
                           V = Vhat,
                           alpha = 1)

p <- mash_compute_posterior_matrices(m$fitted_g, 
                                     mash.data, 
                                     algorithm.version = "Rcpp")

saveRDS(p,file = paste0("/Volumes/mesbah-1/projects/Meta_GWAS/MetaGWAS_N900k/mash/chip_mash_batch",i,".oct2024.rds"))

##########  
file <- readRDS("/Volumes/mesbah-1/projects/Meta_GWAS/MetaGWAS_N900k/mash/chip_mash_batch0.oct2024.rds")

pm <- file$PosteriorMean

lf <- file$lfsr

se <- file$PosteriorSD

for(i in 1:2){
  
  print(i)
  
  file <- readRDS(paste0("/Volumes/mesbah-1/projects/Meta_GWAS/MetaGWAS_N900k/mash/chip_mash_batch",i,".oct2024.rds"))
  
  p <- file$PosteriorMean
  
  l <- file$lfsr
  
  s <- file$PosteriorSD
  
  pm <- rbind(pm, p)
  
  lf <- rbind(lf, l)
  
  se <- rbind(se,s)

}

saveRDS(se,
        "/Volumes/mesbah-1/projects/Meta_GWAS/MetaGWAS_N900k/mash/posterior_sd.oct2024.rds")

saveRDS(pm,
        "/Volumes/mesbah-1/projects/Meta_GWAS/MetaGWAS_N900k/mash/posterior_means.oct2024.rds")

saveRDS(lf,
        "/Volumes/mesbah-1/projects/Meta_GWAS/MetaGWAS_N900k/mash/posterior_lfsr.oct2024.rds")

ptab <- 2*pnorm(-abs(zmash))

sum(lf<0.05)

### 
save.image("Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/chip3_mpn_ltl_mca.mash.oct2024.rda")
###

###### Data driven method of MASH
# https://stephenslab.github.io/mashr/articles/intro_mash_dd.html

mash_plot_meta(m,1)
mash_plot_meta(m,2)
mash_plot_meta(m,1)
mash_plot_meta(m,1)
mash_plot_meta(m,1)


## Plot ### 
rm(list = ls())
load("~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/chip3_mpn_ltl_mca.mash.oct2024.rda")

library("lattice")

pdf("../../../mash/mash.corMtx.pdf", 
    width = 8.5, height = 8)
for(i in 1:7){
  z.num <- as.matrix(cov2cor(m$fitted_g$Ulist[[i]]))
  colnames(z.num) <- colnames(zmash)[c(1:6)]
  rownames(z.num) <- NULL
  clrs = colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#FFFFBF", "#E0F3F8","#91BFDB","#4575B4")))(64)
  z.num[lower.tri(z.num)] <- NA
  print(levelplot(z.num,col.regions = clrs,xlab = "",ylab = "",colorkey = TRUE,main=paste0(names(m$fitted_g$Ulist)[[i]])))
}
dev.off()

## Posterior mean
head(pm)
head(lf)
sig_lfsr <- rowSums(lf<=0.05)
hist(sig_lfsr[sig_lfsr>1], freq = FALSE, main = "Number of Conditions")
## 
library(data.table)
posterior_mean <- as.data.table(pm)
posterior_mean$VarID <- row.names(pm)
head(posterior_mean)
names(posterior_mean) <- c("CHIP PM", "DNMT3A PM", 
                           "TET2 PM","LTL PM", 
                           "mCA PM", "MPN PM", "varID")

# LFSR
LFSR <- as.data.table(lf)
LFSR$VarID <- row.names(lf) 
head(LFSR)
names(LFSR) <- c("CHIP LFSR", "DNMT3A LFSR", 
                 "TET2 LFSR","LTL LFSR", 
                 "mCA LFSR", "MPN LFSR", "varID")
# GWAS P
gwas_p <- as.data.table(ptab)
gwas_p$varID <- row.names(ptab) 
head(gwas_p)
names(gwas_p) <- c("CHIP GWAS P", "DNMT3A GWAS P", 
                 "TET2 GWAS P","LTL GWAS P", 
                 "mCA GWAS P", "MPN GWAS P", "varID")
# combine mash and GWAS result
res_mash <- merge(posterior_mean, LFSR, by="varID")
head(res_mash)
res_mash <- merge(res_mash, gwas_p, by="varID")
head(res_mash)
names(res_mash)
# CHIP
cat("CHIP mash=",sum(res_mash[[8]]<=0.05),
    "GWAS=",sum(res_mash[[14]]<=5e-8),
    ", Fold change=",sum(res_mash[[8]]<=0.05)/sum(res_mash[[14]]<=5e-8),"\n")
# CHIP mash= 21519 GWAS= 2121 , Fold change= 10.14569 
# DNMT3A
cat("DNMT3A mash=",sum(res_mash[[9]]<=0.05),
    "GWAS=",sum(res_mash[[15]]<=5e-8),
    ", Fold change=",sum(res_mash[[9]]<=0.05)/sum(res_mash[[15]]<=5e-8),"\n")
# DNMT3A mash= 20312 GWAS= 2654 , Fold change= 7.653353  
# TET2
cat("TET2 mash=",sum(res_mash[[10]]<=0.05),
    "GWAS=",sum(res_mash[[16]]<=5e-8),
    ", Fold change=",sum(res_mash[[10]]<=0.05)/sum(res_mash[[16]]<=5e-8),"\n")
# TET2 mash= 16169 GWAS= 356 , Fold change= 45.41854 
# LTL
cat("LTL mash=",sum(res_mash[[11]]<=0.05),
    "GWAS=",sum(res_mash[[17]]<=5e-8),
    ", Fold change=",sum(res_mash[[11]]<=0.05)/sum(res_mash[[17]]<=5e-8),"\n")
# LTL mash= 68147 GWAS= 23650 , Fold change= 2.88148
# mCA
cat("mCA mash=",sum(res_mash[[12]]<=0.05),
    "GWAS=",sum(res_mash[[18]]<=5e-8),
    ", Fold change=",sum(res_mash[[12]]<=0.05)/sum(res_mash[[18]]<=5e-8),"\n")
# mCA mash= 14040 GWAS= 3532 , Fold change= 3.975085 
# MPN
cat("MPN mash=",sum(res_mash[[13]]<=0.05),
    "GWAS=",sum(res_mash[[19]]<=5e-8),
    ", Fold change=",sum(res_mash[[13]]<=0.05)/sum(res_mash[[19]]<=5e-8),"\n")
# MPN mash= 13308 GWAS= 749 , Fold change= 17.76769

#
library(dplyr)

#### lfsr<=0.05 | GWAS P <= 5e-8
res_mash_sig <- res_mash %>% 
  filter(`CHIP LFSR`<=0.05 | 
           `DNMT3A LFSR`<=0.05 | 
           `TET2 LFSR`<=0.05 | 
           `LTL LFSR`<=0.05 | 
           `mCA LFSR`<=0.05 | 
           `MPN LFSR`<=0.05 | 
           `CHIP GWAS P`<=5e-8 | 
           `DNMT3A GWAS P`<=5e-8 | 
           `TET2 GWAS P`<=5e-8 |
           `LTL GWAS P`<=5e-8 | 
           `mCA GWAS P`<=5e-8 | 
           `MPN GWAS P`<=5e-8) 

names(res_mash_sig)
str(res_mash_sig)
names(res_mash_sig[,c(1,2,8,14, 
                      3,9,15, 
                      4,10,16, 
                      5,11,17, 
                      6,12,18,
                      7,13,19)])
###
fwrite(res_mash_sig[,c(1,2,8,14, 
                       3,9,15, 
                       4,10,16, 
                       5,11,17, 
                       6,12,18,
                       7,13,19)], "~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/chip_dnmt_tet_ltl_mCA_mpn.lfsr05.oct2024.csv.gz", 
       row.names = F, sep=",")
## E-10
res_mash_sig$`CHIP PM`<- formatC(x = round(res_mash_sig$`CHIP PM`,3), digits = 3, format = "f")
res_mash_sig$`DNMT3A PM`<- formatC(x = round(res_mash_sig$`DNMT3A PM`,3), digits = 3, format = "f")
res_mash_sig$`TET2 PM`<- formatC(x = round(res_mash_sig$`TET2 PM`,3), digits = 3, format = "f")
res_mash_sig$`LTL PM`<- formatC(x = round(res_mash_sig$`LTL PM`,3), digits = 3, format = "f")
res_mash_sig$`mCA PM`<- formatC(x = round(res_mash_sig$`mCA PM`,3), digits = 3, format = "f")
res_mash_sig$`MPN PM`<- formatC(x = round(res_mash_sig$`MPN PM`,3), digits = 3, format = "f")
str(res_mash_sig)

res_mash_sig$`CHIP LFSR` <- formatC(x = res_mash_sig$`CHIP LFSR`, digits = 1, format = "E")
res_mash_sig$`DNMT3A LFSR` <- formatC(x = res_mash_sig$`DNMT3A LFSR`, digits = 1, format = "E")
res_mash_sig$`TET2 LFSR` <- formatC(x = res_mash_sig$`TET2 LFSR`, digits = 1, format = "E")
res_mash_sig$`LTL LFSR` <- formatC(x = res_mash_sig$`LTL LFSR`, digits = 1, format = "E")
res_mash_sig$`mCA LFSR` <- formatC(x = res_mash_sig$`mCA LFSR`, digits = 1, format = "E")
res_mash_sig$`MPN LFSR` <- formatC(x = res_mash_sig$`MPN LFSR`, digits = 1, format = "E")
str(res_mash_sig)

res_mash_sig$`CHIP GWAS P` <- formatC(x = res_mash_sig$`CHIP GWAS P`, digits = 1, format = "E")
res_mash_sig$`DNMT3A GWAS P` <- formatC(x = res_mash_sig$`DNMT3A GWAS P`, digits = 1, format = "E")
res_mash_sig$`TET2 GWAS P` <- formatC(x = res_mash_sig$`TET2 GWAS P`, digits = 1, format = "E")
res_mash_sig$`LTL GWAS P` <- formatC(x = res_mash_sig$`LTL GWAS P`, digits = 1, format = "E")
res_mash_sig$`mCA GWAS P` <- formatC(x = res_mash_sig$`mCA GWAS P`, digits = 1, format = "E")
res_mash_sig$`MPN GWAS P` <- formatC(x = res_mash_sig$`MPN GWAS P`, digits = 1, format = "E")
str(res_mash_sig)
head(res_mash_sig)

fwrite(res_mash_sig[,c(1,2,8,14, 
                       3,9,15, 
                       4,10,16, 
                       5,11,17, 
                       6,12,18,
                       7,13,19)], "~/Documents/Project/Co-Authors/On-going/CHIP ~ GWAS - Mesbah/MetaGWAS_N900k/mash/chip_dnmt_tet_ltl_mCA_mpn.lfsr05.oct2024.csv", row.names = F, sep=",")

# sorted
# head -1 chip_dnmt_tet_ltl_mCA_mpn.lfsr05.oct2024.csv > sorted.chip_dnmt_tet_ltl_mCA_mpn.lfsr05.oct2024.csv
# sed -e '1d' -e 's:,:\t:g' chip_dnmt_tet_ltl_mCA_mpn.lfsr05.oct2024.csv | sort -V -k1 | sed 's:\t:,:g' >> sorted.chip_dnmt_tet_ltl_mCA_mpn.lfsr05.oct2024.csv 
##
# LFSR$VarID <- gsub(pattern = "X", replacement = "", x = LFSR$varID)
# library(stringr)
# head(LFSR$VarID)
# LFSR$VarID <- paste(str_split_fixed(string = LFSR$VarID, 
#                            pattern = "\\.", n = 4)[,1] , 
#                     str_split_fixed(string = LFSR$VarID,
#                                     pattern = "\\.", n = 4)[,2],
#                     str_split_fixed(string = LFSR$VarID,
#                                     pattern = "\\.", n = 4)[,3],
#                     str_split_fixed(string = LFSR$VarID,
#                                     pattern = "\\.", n = 4)[,4],
#                     sep = ":")
# 
# head(LFSR$VarID)

# change_format <- function(variant){
#   variant <- gsub("^X", "", variant)
#   # Split the string by dots
#   parts <- unlist(strsplit(variant, "\\."))
#   paste(parts[1], parts[2], parts[3], parts[4], sep=":")
# }


CHIPs_lfsr <- subset(LFSR[,c(7,1:3)], 
                     LFSR$overallCH_Z<0.05 | 
                       LFSR$DNMT3A_CH_Z<0.05 | 
                       LFSR$TET2_CH_Z<0.05)

CHIPs_lfsr$VarID <- sapply(CHIPs_lfsr$varID, change_format )

head(CHIPs_lfsr$VarID)

fwrite(CHIPs_lfsr[,c(5,2:4)], "../../../mash/chip_dnmt_tet.lfsr05.csv", row.names = F)

## 
library(qqman)
library(stringr)
chip <- CHIPs_lfsr[,c(5,2)]
chip$chr <- as.numeric(str_split_fixed(string = chip$VarID, 
                                      pattern = ":", n = 4)[,1] )
chip$bp <- as.numeric(str_split_fixed(string = chip$VarID, 
                                      pattern = ":", n = 4)[,2] )
head(chip)
# change 0 values to some small values
chip$overallCH_Z[chip$overallCH_Z==0] <- 5e-200
manhattan(x=chip,chr = "chr", 
          bp = "bp", p = "overallCH_Z", 
          snp = "VarID",
          suggestiveline = FALSE, 
          genomewideline = .05)





