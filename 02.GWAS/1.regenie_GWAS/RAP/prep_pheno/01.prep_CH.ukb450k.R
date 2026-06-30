
library(data.table)

## all 450k WES sample ID
ukb450_eid <- fread("/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/450k/ukb450k.eid.list", 
                    header = F)

## UKB 200k CHIP call: Natarajan Lab 
# 200,629
all_ukb200_sam <- fread("/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb200k_eid.csv.gz", header = T)

  # CH variants
ch_var_ukb200k <- fread("/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/200k/UKB200k_CHIP_Variants.11Aug2021.csv", 
                        header = T)
ch_var_ukb200k$chr_pos <- paste(ch_var_ukb200k$CHR, 
                                ch_var_ukb200k$POS_hg38, 
                                sep = "_")

ch_var_ukb200k$Sample_chr_pos <- paste(ch_var_ukb200k$eid_7089, 
                                       ch_var_ukb200k$chr_pos, 
                                       sep = "_")
  # CHIP status
ch_call_ukb200k <- fread("/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/200k/UKB200k_CHIP_Call_11Aug2021.csv", 
                         header = T)
  # 500 samples used in PON
pon_ukb200k <- all_ukb200_sam$eid_7089[!(all_ukb200_sam$eid_7089 %in% ch_call_ukb200k$eid_7089)]

table(ukb450_eid$V1 %in% all_ukb200_sam$eid_7089)
# FALSE   TRUE 
# 254225 200567


### Bick Lab CHIP call
ab_ukb450k_latest <- fread("/Users/muddin/Documents/Project/CHIP_annotation/2022_CHIP_Call/UKB450k/2022Oct/CHIP_calls_Oct16_2022.txt")
boxplot(ab_ukb450k_latest$AF[ab_ukb450k_latest$minAD>=5], ab_ukb450k_latest$AF[ab_ukb450k_latest$minAD<5])
ab_ukb450k_latest$chr_pos <- paste(ab_ukb450k_latest$Chr, ab_ukb450k_latest$Start, sep="_")
hist(ab_ukb450k_latest$AF)
summary(ab_ukb450k_latest$AF)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0200  0.0530  0.0800  0.1165  0.1400  0.9780 
summary(ab_ukb450k_latest$minAD)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.000   3.000   4.000   7.313   8.000 104.000
summary(ab_ukb450k_latest$AF[ab_ukb450k_latest$minAD>=5]) 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0200  0.0870  0.1330  0.1691  0.2170  0.9780 
summary(ab_ukb450k_latest$AF[ab_ukb450k_latest$minAD<5])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.02000 0.04300 0.05500 0.06258 0.07300 0.27500
  # ukb200k subset
ab_ukb200k <- subset(ab_ukb450k_latest, 
                     ab_ukb450k_latest$Broad_ID %in% all_ukb200_sam$eid_7089)


# ukb250k
ab_ukb250k <- subset(ab_ukb450k_latest, 
                     !(ab_ukb450k_latest$Broad_ID %in% all_ukb200_sam$eid_7089) )


## filter U2AF1 
ab_ukb200k_u2af1 <- subset(ab_ukb200k, 
                           ab_ukb200k$Gene.refGene=="U2AF1" & 
                             ab_ukb200k$minAD>=5)

ab_ukb250k_u2af1 <- subset(ab_ukb250k, 
                           ab_ukb250k$Gene.refGene=="U2AF1" & 
                             ab_ukb250k$minAD>=5)
ab_ukb250k_ETNK1 <- subset(ab_ukb250k, 
                           ab_ukb250k$Gene.refGene=="ETNK1")

ab_ukb250k_ETNK1 <- ab_ukb250k_ETNK1[,c(1,3,4,7,8,34,31,32,33,2,28,29)]
names(ab_ukb250k_ETNK1) <- names(ch_var_ukb200k_pls_u2af1)

ab_ukb250k_no_u2af1_ETNK1 <- subset(ab_ukb250k, 
                                    !(ab_ukb250k$Gene.refGene %in% c("U2AF1", 
                                                                     "ETNK1") ))
ab_ukb250k_no_u2af1_ETNK1$F1R2_2 <- as.numeric(ab_ukb250k_no_u2af1_ETNK1$F1R2_2 )
ab_ukb250k_no_u2af1_ETNK1$F2R1_2 <- as.numeric(ab_ukb250k_no_u2af1_ETNK1$F2R1_2 )
table(ab_ukb250k_no_u2af1_ETNK1$F1R2_2>=1 & 
        ab_ukb250k_no_u2af1_ETNK1$F2R1_2>=1)
# FALSE  TRUE 
# 14 19266
# FR/RR>=1
ab_ukb250k_no_u2af1_ETNK1 <- ab_ukb250k_no_u2af1_ETNK1[(ab_ukb250k_no_u2af1_ETNK1$F1R2_2>=1 & 
                                                          ab_ukb250k_no_u2af1_ETNK1$F2R1_2>=1),c(1,16,17,19,20,35,31,32,33,2,28,29)]
names(ab_ukb250k_no_u2af1_ETNK1) <- names(ab_ukb250k_ETNK1)

ab_ukb250k_filtered <- as.data.frame(rbind(ab_ukb250k_no_u2af1_ETNK1,
                                           ab_ukb250k_ETNK1,
                                           ab_ukb250k_u2af1))

table(table(ab_ukb250k_filtered$eid_7089))
# 1     2     3     4     5     6 
# 16633  1118   113    18     2     3
ab_ukb250k_filtered$Batch <- "ukb250k"

## PN lab + U2AF1 from Bick Lab
names(ch_var_ukb200k)
ch_var_ukb200k_pls_u2af1 <- ch_var_ukb200k[,c(2,5:12,18,21,28)]

names(ab_ukb200k_u2af1)
ab_ukb200k_u2af1 <- ab_ukb200k_u2af1[,c(1,3,4,6,7,35,31,32,33,2,11,29)]
names(ab_ukb200k_u2af1) <- names(ch_var_ukb200k_pls_u2af1)

ab_ukb250k_u2af1 <- ab_ukb250k_u2af1[,c(1,3,4,6,7,35,31,32,33,2,11,29)]
names(ab_ukb250k_u2af1) <- names(ch_var_ukb200k_pls_u2af1)

ch_var_ukb200k_pls_u2af1 <- as.data.frame(rbind(ch_var_ukb200k_pls_u2af1, ab_ukb200k_u2af1))
ch_var_ukb200k_pls_u2af1$Batch <- "ukb200k"

# remove "p."
ch_var_ukb200k_pls_u2af1$Protein_Change <- gsub(pattern = "p.", 
                                                replacement = "",
                                                x = ch_var_ukb200k_pls_u2af1$Protein_Change )

## all ch variants in ukb200k and ukb250k
ch_var_all_ukb200k_n_ukb250k <- as.data.frame(rbind(ch_var_ukb200k_pls_u2af1, 
                                                    ab_ukb250k_filtered))
head(prop.table(sort(table(ch_var_all_ukb200k_n_ukb250k$Gene),decreasing = T))*100, 10)
# DNMT3A       TET2      ASXL1      PPM1D       TP53      SRSF2 
# 53.0328325 18.4597786  8.6687689  2.8720707  2.2444815  1.1407902 
# JAK2      SF3B1      ASXL2     ZNF318 
# 1.0480430  0.9058307  0.6894206  0.5904903 
summary(ch_var_all_ukb200k_n_ukb250k$AF)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0100  0.0590  0.0910  0.1312  0.1610  0.9780
# VAF>=.02
ch_var_all_ukb200k_n_ukb250k <- subset(ch_var_all_ukb200k_n_ukb250k, 
                                       ch_var_all_ukb200k_n_ukb250k$AF>=0.02)

## all ukb CHIP VAF>=2%
# fwrite(ch_var_all_ukb200k_n_ukb250k, "/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/ch_var_all_ukb200k_n_ukb250k.7Mar2023.csv.gz", 
#        row.names = F, col.names = T, sep=",", na = "NA")

## CH status "0/1"
ch_var_all_ukb200k_n_ukb250k <- fread("/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/ch_var_all_ukb200k_n_ukb250k.7Mar2023.csv.gz", 
                                      header = T)

names(ukb450_eid) <- "eid_7089"

ukb450_eid$Batch <- ifelse(ukb450_eid$eid_7089 %in% all_ukb200_sam$eid_7089, "ukb200k","ukb250k")
ukb450_eid$PON <- ifelse(ukb450_eid$eid_7089 %in% pon_ukb200k, 1, 0)
# VAF2
ukb450_eid$hasCH <- ifelse(ukb450_eid$eid_7089 %in% unique(ch_var_all_ukb200k_n_ukb250k$eid_7089), 1, 0)
prop.table(table(ukb450_eid$hasCH))
# 0          1 
# 0.93475479 0.06524521 
# VAF 5
ukb450_eid$hasCHvaf05 <- ifelse(ukb450_eid$eid_7089 %in% unique(ch_var_all_ukb200k_n_ukb250k$eid_7089[ch_var_all_ukb200k_n_ukb250k$AF>=0.05]), 1, 0)
prop.table(table(ukb450_eid$hasCHvaf05))
# 0          1 
# 0.94399198 0.05600802
# VAF 10
ukb450_eid$hasCHvaf10 <- ifelse(ukb450_eid$eid_7089 %in% unique(ch_var_all_ukb200k_n_ukb250k$eid_7089[ch_var_all_ukb200k_n_ukb250k$AF>=0.10]), 1, 0)
prop.table(table(ukb450_eid$hasCHvaf10))
# 0          1 
# 0.97029191 0.02970809 
table(ukb450_eid$hasCHvaf10, ukb450_eid$hasCH, exclude = NULL)

ukb450_eid$hasDNMT3A <- ifelse(ukb450_eid$eid_7089 %in% unique(ch_var_all_ukb200k_n_ukb250k$eid_7089[ch_var_all_ukb200k_n_ukb250k$Gene=="DNMT3A"]), 1, ifelse(ukb450_eid$hasCH==0,0,NA))

table(ukb450_eid$hasDNMT3A, ukb450_eid$hasCH, exclude = NULL)

table(ukb450_eid$hasDNMT3A, exclude = NULL)
# 0      1   <NA> 
#   425119  16605  13068 

ukb450_eid$hasTET2 <- ifelse(ukb450_eid$eid_7089 %in% unique(ch_var_all_ukb200k_n_ukb250k$eid_7089[ch_var_all_ukb200k_n_ukb250k$Gene=="TET2"]), 1, ifelse(ukb450_eid$hasCH==0,0,NA))

table(ukb450_eid$hasTET2, exclude = NULL)
# 0      1   <NA> 
#   425119   5710  23963
ukb450_eid$hasASXL1 <- ifelse(ukb450_eid$eid_7089 %in% unique(ch_var_all_ukb200k_n_ukb250k$eid_7089[ch_var_all_ukb200k_n_ukb250k$Gene=="ASXL1"]), 1, ifelse(ukb450_eid$hasCH==0,0,NA))

table(ukb450_eid$hasASXL1, exclude = NULL)
# 0      1   <NA> 
#   425119   2775  26898
## ASXL1, DNMT3A , TET2
ukb450_eid$hasDTA <- ifelse(ukb450_eid$eid_7089 %in% unique(ch_var_all_ukb200k_n_ukb250k$eid_7089[ch_var_all_ukb200k_n_ukb250k$Gene=="ASXL1" | ch_var_all_ukb200k_n_ukb250k$Gene== "DNMT3A" | ch_var_all_ukb200k_n_ukb250k$Gene== "TET2"]), 1, ifelse(ukb450_eid$hasCH==0,0,NA))

table(ukb450_eid$hasDTA, exclude = NULL)
# 0      1   <NA> 
#   425119  24289   5384
## Splicing factors: SF3B1, U2AF1, SRSF2, ZRSR2
#c("SRSF2","SF3B1","U2AF1","ZRSR2")
ukb450_eid$hasSF <- ifelse(ukb450_eid$eid_7089 %in% unique(ch_var_all_ukb200k_n_ukb250k$eid_7089[ch_var_all_ukb200k_n_ukb250k$Gene %in% c("SRSF2","SF3B1","U2AF1","ZRSR2")] ), 1, ifelse(ukb450_eid$hasCH==0,0,NA))

table(ukb450_eid$hasSF, exclude = NULL)
# 0      1   <NA> 
#   425119    736  28937

## DDR: TP53   PPM1D
ukb450_eid$hasDDR <- ifelse(ukb450_eid$eid_7089 %in% unique(ch_var_all_ukb200k_n_ukb250k$eid_7089[ch_var_all_ukb200k_n_ukb250k$Gene %in% c("PPM1D","TP53")] ), 1, ifelse(ukb450_eid$hasCH==0,0,NA))

table(ukb450_eid$hasDDR, exclude = NULL)
# 0      1   <NA> 
#   425119   1632  28041

# fwrite(ukb450_eid, "/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/ch_call_ukb200k_n_ukb250k.7Mar2023.csv.gz",
#        row.names = F, col.names = T, sep=",", na = "NA")

# save.image(file = "/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/curation_ukb200k_n_ukb250k.7Mar2023.rda")
