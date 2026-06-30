###########################################################
################## ABIC Lab (Vlasschaert, Caitlyn ) ############
library(readxl)
library(data.table)
library(dplyr)
library(tidyr)
var_chip.ukb450 <- read_xlsx("~/Documents/Project/CHIP_annotation/2022_CHIP_Call/UKB450k/CHIP_call_set_07182022.xlsx", 
                             sheet = 2)

ukb450k.v1 <- var_chip.ukb450 %>% filter(Otherinfo12== "GT:AD:AF:DP:F1R2:F2R1:SB") %>% separate(Otherinfo13, c("GT","AD","AF","DP","F1R2","F2R1","SB"), ":")
ukb450k.v2 <- var_chip.ukb450 %>% filter(Otherinfo12=="GT:AD:AF:DP:F1R2:F2R1:PGT:PID:PS:SB") %>% separate(Otherinfo13, c("GT","AD","AF","DP","F1R2","F2R1","PGT","PID","PS","SB"), ":")
gc()
ukb450k.v2$PGT=NULL
ukb450k.v2$PID=NULL
ukb450k.v2$PS=NULL
ukb450k.v2$SB=NULL
ukb450k.v1$SB=NULL
ukb450k <- as.data.frame(rbind(ukb450k.v1, ukb450k.v2)); rm(ukb450k.v1, ukb450k.v2)
ukb450k$FR.Ref <- as.numeric( stringr::str_split_fixed(string = ukb450k$F1R2, pattern = "[,]", n = 2)[,1])
ukb450k$FR.Alt <- as.numeric( stringr::str_split_fixed(string = ukb450k$F1R2, pattern = "[,]", n = 2)[,2])
ukb450k$RR.Ref <- as.numeric( stringr::str_split_fixed(string = ukb450k$F2R1, pattern = "[,]", n = 2)[,1])
ukb450k$RR.Alt <- as.numeric( stringr::str_split_fixed(string = ukb450k$F2R1, pattern = "[,]", n = 2)[,2])
ukb450k$AD.Alt <- as.numeric( stringr::str_split_fixed(string = ukb450k$AD, pattern = "[,]", n = 2)[,2])
ukb450k$DP <- as.numeric(ukb450k$DP)
ukb450k$VAF <- as.numeric(ukb450k$AF)
# ukb450k$gnomadAF <- as.numeric(ukb450k$AF_raw)
ukb450k$varID <- paste(ukb450k$Otherinfo4, ukb450k$Otherinfo5, 
                       ukb450k$Otherinfo7, ukb450k$Otherinfo8, 
                       sep = "_")

summary(ukb450k$AD.Alt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.000   3.000   4.000   7.301   8.000 104.000
summary(ukb450k$FR.Alt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   1.000   2.000   3.606   4.000  48.000
summary(ukb450k$RR.Alt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    1.00    2.00    3.55    4.00   55.00
ukb450k <- subset(ukb450k, ukb450k$FR.Alt>=1)
cor(ukb450k$VAF, ukb450k$AD.Alt/ukb450k$DP)
# 0.9899306
plot(ukb450k$VAF, ukb450k$AD.Alt/ukb450k$DP)
## Exclude non-CHIP genes
# not_chip_gene <- c("SF3A1","CSF3R", "PRPF40B", "SF1", "PDSS2")

var_chip.ukb450.CHIPonly <- subset(ukb450k, !(ukb450k$Gene.refGene %in% c("SF3A1","CSF3R", "PRPF40B", "SF1", "PDSS2")) )

summary(var_chip.ukb450.CHIPonly)

# fwrite(var_chip.ukb450.CHIPonly, "/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/chipVar_CV_AB.ukb450k_all.csv.gz",
#        row.names = F, col.names = T, sep=",", na = 'NA')

