
################################ Variant Allele Fraction #########
################################ 
library(data.table)
############ MGBB 53k ################
## variant allele fractions
mgbb_ch <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/chip_call/mgbb53k.chip_call_all52966.tsv.gz")
mgbb_ch <- mgbb_ch[,c(1:6,27,7:10,14,15,28:47)]

mgbb_var <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.chip_var_N5910_nchip7145.csv.gz") 
sort(table(mgbb_var$Gene), decreasing = T)

# Annotate 
mgbb_ch$CHIP <- ifelse(mgbb_ch$Biobank_Subject_ID %in% mgbb_var$sample_id[mgbb_var$VAF>=0.02], 1, 0)
mgbb_ch$DNMT3A <- ifelse(mgbb_ch$Biobank_Subject_ID %in% mgbb_var$sample_id[mgbb_var$VAF>=0.02 & mgbb_var$Gene=="DNMT3A"], 1, ifelse(mgbb_ch$CHIP==0, 0, NA))
mgbb_ch$TET2 <- ifelse(mgbb_ch$Biobank_Subject_ID %in% mgbb_var$sample_id[mgbb_var$VAF>=0.02 & mgbb_var$Gene=="TET2"], 1, ifelse(mgbb_ch$CHIP==0, 0, NA))
mgbb_ch$ASXL1 <- ifelse(mgbb_ch$Biobank_Subject_ID %in% mgbb_var$sample_id[mgbb_var$VAF>=0.02 & mgbb_var$Gene=="ASXL1"], 1, ifelse(mgbb_ch$CHIP==0, 0, NA))
mgbb_ch$JAK2 <- ifelse(mgbb_ch$Biobank_Subject_ID %in% mgbb_var$sample_id[mgbb_var$VAF>=0.02 & mgbb_var$Gene=="JAK2"], 1, ifelse(mgbb_ch$CHIP==0, 0, NA))
mgbb_ch$PPM1D <- ifelse(mgbb_ch$Biobank_Subject_ID %in% mgbb_var$sample_id[mgbb_var$VAF>=0.02 & mgbb_var$Gene=="PPM1D"], 1, ifelse(mgbb_ch$CHIP==0, 0, NA))
mgbb_ch$TP53 <- ifelse(mgbb_ch$Biobank_Subject_ID %in% mgbb_var$sample_id[mgbb_var$VAF>=0.02 & mgbb_var$Gene=="TP53"], 1, ifelse(mgbb_ch$CHIP==0, 0, NA))
mgbb_ch$SF3B1 <- ifelse(mgbb_ch$Biobank_Subject_ID %in% mgbb_var$sample_id[mgbb_var$VAF>=0.02 & mgbb_var$Gene=="SF3B1"], 1, ifelse(mgbb_ch$CHIP==0, 0, NA))
mgbb_ch$SRSF2 <- ifelse(mgbb_ch$Biobank_Subject_ID %in% mgbb_var$sample_id[mgbb_var$VAF>=0.02 & mgbb_var$Gene=="SRSF2"], 1, ifelse(mgbb_ch$CHIP==0, 0, NA))
mgbb_ch$DTA <- ifelse(mgbb_ch$Biobank_Subject_ID %in% mgbb_var$sample_id[mgbb_var$VAF>=0.02 & (mgbb_var$Gene %in% c("DNMT3A","TET2","ASXL1"))], 1, ifelse(mgbb_ch$CHIP==0, 0, NA))
mgbb_ch$DDR <- ifelse(mgbb_ch$Biobank_Subject_ID %in% mgbb_var$sample_id[mgbb_var$VAF>=0.02 & (mgbb_var$Gene %in% c("PPM1D","TP53"))], 1, ifelse(mgbb_ch$CHIP==0, 0, NA))
mgbb_ch$Splice <- ifelse(mgbb_ch$Biobank_Subject_ID %in% mgbb_var$sample_id[mgbb_var$VAF>=0.02 & (mgbb_var$Gene %in% c("SF3B1","SRSF2"))], 1, ifelse(mgbb_ch$CHIP==0, 0, NA))

mgbb_ch$expandedCHIP <- ifelse(mgbb_ch$Biobank_Subject_ID %in% mgbb_var$sample_id[mgbb_var$VAF>=0.10], 1, ifelse(mgbb_ch$CHIP==0, 0, NA))
mgbb_ch$expandedDNMT3A <- ifelse(mgbb_ch$Biobank_Subject_ID %in% mgbb_var$sample_id[mgbb_var$VAF>=0.10 & mgbb_var$Gene=="DNMT3A"], 1, ifelse(mgbb_ch$CHIP==0, 0, NA))
mgbb_ch$expandedTET2 <- ifelse(mgbb_ch$Biobank_Subject_ID %in% mgbb_var$sample_id[mgbb_var$VAF>=0.1 & mgbb_var$Gene=="TET2"], 1, ifelse(mgbb_ch$CHIP==0, 0, NA))
mgbb_ch$expandedASXL1 <- ifelse(mgbb_ch$Biobank_Subject_ID %in% mgbb_var$sample_id[mgbb_var$VAF>=0.1 & mgbb_var$Gene=="ASXL1"], 1, ifelse(mgbb_ch$CHIP==0, 0, NA))
mgbb_ch$expandedJAK2 <- ifelse(mgbb_ch$Biobank_Subject_ID %in% mgbb_var$sample_id[mgbb_var$VAF>=0.1 & mgbb_var$Gene=="JAK2"], 1, ifelse(mgbb_ch$CHIP==0, 0, NA))
mgbb_ch$expandedPPM1D <- ifelse(mgbb_ch$Biobank_Subject_ID %in% mgbb_var$sample_id[mgbb_var$VAF>=0.1 & mgbb_var$Gene=="PPM1D"], 1, ifelse(mgbb_ch$CHIP==0, 0, NA))
mgbb_ch$expandedTP53 <- ifelse(mgbb_ch$Biobank_Subject_ID %in% mgbb_var$sample_id[mgbb_var$VAF>=0.1 & mgbb_var$Gene=="TP53"], 1, ifelse(mgbb_ch$CHIP==0, 0, NA))
mgbb_ch$expandedSF3B1 <- ifelse(mgbb_ch$Biobank_Subject_ID %in% mgbb_var$sample_id[mgbb_var$VAF>=0.1 & mgbb_var$Gene=="SF3B1"], 1, ifelse(mgbb_ch$CHIP==0, 0, NA))
mgbb_ch$expandedSRSF2 <- ifelse(mgbb_ch$Biobank_Subject_ID %in% mgbb_var$sample_id[mgbb_var$VAF>=0.1 & mgbb_var$Gene=="SRSF2"], 1, ifelse(mgbb_ch$CHIP==0, 0, NA))
mgbb_ch$expandedDTA <- ifelse(mgbb_ch$Biobank_Subject_ID %in% mgbb_var$sample_id[mgbb_var$VAF>=0.1 & (mgbb_var$Gene %in% c("DNMT3A","TET2","ASXL1"))], 1, ifelse(mgbb_ch$CHIP==0, 0, NA))
mgbb_ch$expandedDDR <- ifelse(mgbb_ch$Biobank_Subject_ID %in% mgbb_var$sample_id[mgbb_var$VAF>=0.1 & (mgbb_var$Gene %in% c("PPM1D","TP53"))], 1, ifelse(mgbb_ch$CHIP==0, 0, NA))
mgbb_ch$expandedSplice <- ifelse(mgbb_ch$Biobank_Subject_ID %in% mgbb_var$sample_id[mgbb_var$VAF>=0.1 & (mgbb_var$Gene %in% c("SF3B1","SRSF2"))], 1, ifelse(mgbb_ch$CHIP==0, 0, NA))

for(i in 34:57){cat(names(mgbb_ch)[i],prop.table(table(mgbb_ch[[i]], exclude = NULL))*100,"\n")}

######  
# fwrite(mgbb_ch,"/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/chip_call/mgbb53k.CHIP_phewas.csv.gz",
#        row.names=F, col.names=T, sep=",", na="NA")
# 
# fwrite(mgbb_var,"/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/chip_call/mgbb53k.variants_phewas.csv.gz",
#        row.names=F, col.names=T, sep=",", na="NA")
######################################


################ UKB450k #################
## UKBB 200k
# ukb200k_eid <- fread("/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/200k/UKB200k_eid_all_wes_analysed.csv", header = T)
  # N=200,128 samples (excluded 500 used in Panel-of-normal)
ukb200k_ch <- fread("/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/200k/UKB200k_CHIP_Call_11Aug2021.csv", header = T)
ukb200k_ch$CHIP_Batch <- "UKB200k"

# UKBB 250k; N=254,225 samples
ukb250k_ch <- fread("/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/450k/chipCall_CV_AB.ukb250k_all.csv.gz", header = T)
ukb250k_ch$CHIP_Batch <- "UKB250k"

  # merge UKB200k and ukb250k data
#  N = 454,353
ukb450k_ch <- as.data.frame(rbind(ukb200k_ch[,c(1,4,6,13)],
                                  ukb250k_ch[,c(3,6,5,44)],
                                  use.names=FALSE)); rm(ukb200k_ch,ukb250k_ch)
ukb450k_ch$SEX[ukb450k_ch$SEX==0] <- "Female"
ukb450k_ch$SEX[ukb450k_ch$SEX==1] <- "Male"

## Exclude Consent withdrawn
ukb_withdrawn <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/w7089_20220222.csv")
ukb450k_ch <- subset(ukb450k_ch, !(ukb450k_ch$eid_7089 %in% ukb_withdrawn$V1)) ; rm(ukb_withdrawn)
nrow(ukb450k_ch)
# 454,327

  # Variants
ukb200k_var <- fread("/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/200k/UKB200k_CHIP_Variants.11Aug2021.csv", header = T)
ukb200k_var$VAF <- ukb200k_var$AF
ukb200k_var$CHIP_Batch <- "UKB200k"

ukb250k_var <- fread("/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/450k/chipVars_CV_AB.ukb250k_all.csv.gz", header = T) 
ukb250k_var$CHIP_Batch <- "UKB250k"
ukb450k_var <- as.data.frame(rbind(ukb200k_var[,c(2,18,39,40)],
                                   ukb250k_var[,c(1,2,54,56)], 
                                  use.names=FALSE)); rm(ukb200k_var,ukb250k_var)

# Annotate 
ukb450k_ch$CHIP <- ifelse(ukb450k_ch$eid_7089 %in% ukb450k_var$eid_7089[ukb450k_var$VAF>=0.02], 1, 0)
ukb450k_ch$DNMT3A <- ifelse(ukb450k_ch$eid_7089 %in% ukb450k_var$eid_7089[ukb450k_var$VAF>=0.02 & ukb450k_var$Gene=="DNMT3A"], 1, ifelse(ukb450k_ch$CHIP==0, 0, NA))
ukb450k_ch$TET2 <- ifelse(ukb450k_ch$eid_7089 %in% ukb450k_var$eid_7089[ukb450k_var$VAF>=0.02 & ukb450k_var$Gene=="TET2"], 1, ifelse(ukb450k_ch$CHIP==0, 0, NA))
ukb450k_ch$ASXL1 <- ifelse(ukb450k_ch$eid_7089 %in% ukb450k_var$eid_7089[ukb450k_var$VAF>=0.02 & ukb450k_var$Gene=="ASXL1"], 1, ifelse(ukb450k_ch$CHIP==0, 0, NA))
ukb450k_ch$JAK2 <- ifelse(ukb450k_ch$eid_7089 %in% ukb450k_var$eid_7089[ukb450k_var$VAF>=0.02 & ukb450k_var$Gene=="JAK2"], 1, ifelse(ukb450k_ch$CHIP==0, 0, NA))
ukb450k_ch$PPM1D <- ifelse(ukb450k_ch$eid_7089 %in% ukb450k_var$eid_7089[ukb450k_var$VAF>=0.02 & ukb450k_var$Gene=="PPM1D"], 1, ifelse(ukb450k_ch$CHIP==0, 0, NA))
ukb450k_ch$TP53 <- ifelse(ukb450k_ch$eid_7089 %in% ukb450k_var$eid_7089[ukb450k_var$VAF>=0.02 & ukb450k_var$Gene=="TP53"], 1, ifelse(ukb450k_ch$CHIP==0, 0, NA))
ukb450k_ch$SF3B1 <- ifelse(ukb450k_ch$eid_7089 %in% ukb450k_var$eid_7089[ukb450k_var$VAF>=0.02 & ukb450k_var$Gene=="SF3B1"], 1, ifelse(ukb450k_ch$CHIP==0, 0, NA))
ukb450k_ch$SRSF2 <- ifelse(ukb450k_ch$eid_7089 %in% ukb450k_var$eid_7089[ukb450k_var$VAF>=0.02 & ukb450k_var$Gene=="SRSF2"], 1, ifelse(ukb450k_ch$CHIP==0, 0, NA))
ukb450k_ch$DTA <- ifelse(ukb450k_ch$eid_7089 %in% ukb450k_var$eid_7089[ukb450k_var$VAF>=0.02 & (ukb450k_var$Gene %in% c("DNMT3A","TET2","ASXL1"))], 1, ifelse(ukb450k_ch$CHIP==0, 0, NA))
ukb450k_ch$DDR <- ifelse(ukb450k_ch$eid_7089 %in% ukb450k_var$eid_7089[ukb450k_var$VAF>=0.02 & (ukb450k_var$Gene %in% c("PPM1D","TP53"))], 1, ifelse(ukb450k_ch$CHIP==0, 0, NA))
ukb450k_ch$Splice <- ifelse(ukb450k_ch$eid_7089 %in% ukb450k_var$eid_7089[ukb450k_var$VAF>=0.02 & (ukb450k_var$Gene %in% c("SF3B1","SRSF2"))], 1, ifelse(ukb450k_ch$CHIP==0, 0, NA))

ukb450k_ch$expandedCHIP <- ifelse(ukb450k_ch$eid_7089 %in% ukb450k_var$eid_7089[ukb450k_var$VAF>=0.10], 1, ifelse(ukb450k_ch$CHIP==0, 0, NA))
ukb450k_ch$expandedDNMT3A <- ifelse(ukb450k_ch$eid_7089 %in% ukb450k_var$eid_7089[ukb450k_var$VAF>=0.10 & ukb450k_var$Gene=="DNMT3A"], 1, ifelse(ukb450k_ch$CHIP==0, 0, NA))
ukb450k_ch$expandedTET2 <- ifelse(ukb450k_ch$eid_7089 %in% ukb450k_var$eid_7089[ukb450k_var$VAF>=0.1 & ukb450k_var$Gene=="TET2"], 1, ifelse(ukb450k_ch$CHIP==0, 0, NA))
ukb450k_ch$expandedASXL1 <- ifelse(ukb450k_ch$eid_7089 %in% ukb450k_var$eid_7089[ukb450k_var$VAF>=0.1 & ukb450k_var$Gene=="ASXL1"], 1, ifelse(ukb450k_ch$CHIP==0, 0, NA))
ukb450k_ch$expandedJAK2 <- ifelse(ukb450k_ch$eid_7089 %in% ukb450k_var$eid_7089[ukb450k_var$VAF>=0.1 & ukb450k_var$Gene=="JAK2"], 1, ifelse(ukb450k_ch$CHIP==0, 0, NA))
ukb450k_ch$expandedPPM1D <- ifelse(ukb450k_ch$eid_7089 %in% ukb450k_var$eid_7089[ukb450k_var$VAF>=0.1 & ukb450k_var$Gene=="PPM1D"], 1, ifelse(ukb450k_ch$CHIP==0, 0, NA))
ukb450k_ch$expandedTP53 <- ifelse(ukb450k_ch$eid_7089 %in% ukb450k_var$eid_7089[ukb450k_var$VAF>=0.1 & ukb450k_var$Gene=="TP53"], 1, ifelse(ukb450k_ch$CHIP==0, 0, NA))
ukb450k_ch$expandedSF3B1 <- ifelse(ukb450k_ch$eid_7089 %in% ukb450k_var$eid_7089[ukb450k_var$VAF>=0.1 & ukb450k_var$Gene=="SF3B1"], 1, ifelse(ukb450k_ch$CHIP==0, 0, NA))
ukb450k_ch$expandedSRSF2 <- ifelse(ukb450k_ch$eid_7089 %in% ukb450k_var$eid_7089[ukb450k_var$VAF>=0.1 & ukb450k_var$Gene=="SRSF2"], 1, ifelse(ukb450k_ch$CHIP==0, 0, NA))
ukb450k_ch$expandedDTA <- ifelse(ukb450k_ch$eid_7089 %in% ukb450k_var$eid_7089[ukb450k_var$VAF>=0.1 & (ukb450k_var$Gene %in% c("DNMT3A","TET2","ASXL1"))], 1, ifelse(ukb450k_ch$CHIP==0, 0, NA))
ukb450k_ch$expandedDDR <- ifelse(ukb450k_ch$eid_7089 %in% ukb450k_var$eid_7089[ukb450k_var$VAF>=0.1 & (ukb450k_var$Gene %in% c("PPM1D","TP53"))], 1, ifelse(ukb450k_ch$CHIP==0, 0, NA))
ukb450k_ch$expandedSplice <- ifelse(ukb450k_ch$eid_7089 %in% ukb450k_var$eid_7089[ukb450k_var$VAF>=0.1 & (ukb450k_var$Gene %in% c("SF3B1","SRSF2"))], 1, ifelse(ukb450k_ch$CHIP==0, 0, NA))

for(i in 5:28){cat(names(ukb450k_ch)[i],prop.table(table(ukb450k_ch[[i]], exclude = NULL))*100,"\n")}



# fwrite(ukb450k_ch,"/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/chip_call/ukb450k.CHIP_phewas.csv.gz",
#        row.names=F, col.names=T, sep=",", na="NA")
# 
# fwrite(ukb450k_var,"/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/chip_call/ukb450k.variants_phewas.csv.gz",
#        row.names=F, col.names=T, sep=",", na="NA")

###########################################
# save.image(file = "/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/chip_call/ukb450k_mgbb53k.CHIP_vars.for_phewas.rda")
###########################################
load("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/chip_call/ukb450k_mgbb53k.CHIP_vars.for_phewas.rda")
