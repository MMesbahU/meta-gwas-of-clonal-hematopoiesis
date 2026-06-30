

  ## Prepare phenotypes for GWAS 
library(data.table)
##
## CHIP call
ch_var_ukb450 <-  fread("/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/ch_var_all_ukb200k_n_ukb250k.7Mar2023.csv.gz", 
                        header = T)

ch_ukb450 <-  fread("/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/ch_call_ukb200k_n_ukb250k.7Mar2023.csv.gz", 
                    header = T)
##
  ## UKB demographic informations
# Load UKB data
# "d_base"    "d_pcs"     "d_ukb_all"
load("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/ukb500k_covariates.rda")
rm(d_base, d_pcs)



d_ukb_450k <- subset(d_ukb_all, 
                     d_ukb_all$eid %in% ch_ukb450$eid_7089); rm(d_ukb_all)

names(d_ukb_450k) <- c("eid", "GT_batch.p22000", 
                       "Age_at_recruitment.p21022",
                       "Sex.self_rep.p31","Sex.Genetic.p22001",
                       "aneuploidy.p22019", 
                       "Genetic_ethnic.p22006",
                       "self_ethnic.p21000_i0",
                       "ever_smoked.p20160_i0",         
                       "Smoking_status.p20116_i0",
                       "BMI.p21001_i0","BMI.p23104_i0",
                       names(d_ukb_450k)[13:26])

################################### annotate CHIP
head(sort(table(c(ch_var_ukb450$Gene) ), decreasing = T),15)
# DNMT3A   TET2  ASXL1  PPM1D   TP53  SRSF2   JAK2  SF3B1  ASXL2 ZNF318 
# 17154   5970   2803    929    726    369    339    293    222    189 
# GNB1    NF1    CBL CREBBP  YLPM1 
# 179    174    166    158    158

## CHIP
d_ukb_450k <- merge(d_ukb_450k, ch_ukb450, 
                    by.x="eid", by.y="eid_7089")

### 
d_ukb_450k$FID <- d_ukb_450k$eid
d_ukb_450k$IID <- d_ukb_450k$eid
d_ukb_450k <- d_ukb_450k[, c(38,39,1:37)]

# Ancestry: British-Caucasian coded as "WB" (White British) vs. "Others"  
d_ukb_450k$Ethnic_Background <- ifelse(d_ukb_450k$self_ethnic.p21000_i0=="British" & 
                                         d_ukb_450k$Genetic_ethnic.p22006=="Caucasian", 
                                       "WB", "Others")
table(d_ukb_450k$Ethnic_Background)
# Others     WB 
# 73686 381106
# Age^2
d_ukb_450k$sqrtAge_at_recruitment <- d_ukb_450k$Age_at_recruitment^2
str(d_ukb_450k)
## 0 sample withdrawn 
ukb_withdrawn <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/w7089_20220222.csv")
table(d_ukb_450k$eid %in% ukb_withdrawn$V1, exclude = NULL)
# FALSE 
# 454792
#################################################### 

### Step by step sample filtering for GWAS
## UKB450k samples N=454792
## 1. PON 500: 454792 - 500 = 454292
d_ukb_450k.nomiss <- subset(d_ukb_450k, 
                            d_ukb_450k$PON==0)
# 2. no PCs 393: 454292 - 393 = 453899
d_ukb_450k.nomiss <- subset(d_ukb_450k.nomiss, 
                            !is.na(d_ukb_450k.nomiss$PC1) )

# 3. Sex mismatch 94: 453899 - 94 = 453805
d_ukb_450k.nomiss <- subset(d_ukb_450k.nomiss, 
                            !is.na(d_ukb_450k.nomiss$Genetic_Sex))

## 4. Aneuploidy 430: 453805 - 430 = 453375
d_ukb_450k.nomiss <- subset(d_ukb_450k.nomiss, 
                            d_ukb_450k.nomiss$Aneuploidy == "No_Aneuploidy" )
# withdrawn 0: 453375 - 0 = 453375
# table(d_ukb_450k.nomiss$eid %in% ukb_withdrawn$V1)
# FALSE 
# 453375

## Exclude related samples
geno_related <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/CHIP_GWAS/PhenoCov_prep/ukb7089_rel_s488264.dat", 
                      header = T, stringsAsFactors = F)

library(ukbtools)

relSample_filter.200k <- ukb_gen_samples_to_remove(
  data = geno_related,
  ukb_with_data = d_ukb_450k.nomiss$eid[d_ukb_450k.nomiss$Batch=="ukb200k"],
  cutoff = 0.0884)

relSample_filter.250k <- ukb_gen_samples_to_remove(
  data = geno_related,
  ukb_with_data = d_ukb_450k.nomiss$eid[d_ukb_450k.nomiss$Batch=="ukb250k"],
  cutoff = 0.0884)

  # VAF samples
relSample_filter.VAF_samples <- ukb_gen_samples_to_remove(
  data = geno_related,
  ukb_with_data = unique(ch_var_ukb450$eid_7089),
  cutoff = 0.0884)

library(plinkQC)
# 
geno_related_sam <- relatednessFilter(geno_related, 
                                      relatednessTh = 0.0884,
                                      relatednessIID1 = "ID1",
                                      relatednessIID2 = "ID2",
                                      relatednessRelatedness = "Kinship",
                                      verbose = FALSE)
sam2remove <- geno_related_sam$failIDs

table(d_ukb_450k.nomiss$eid %in% sam2remove$IID)
# FALSE   TRUE 
# 421316  32059
  # UKB200k only
d_ukb200k.nomiss <- subset(d_ukb_450k.nomiss, 
                           d_ukb_450k.nomiss$Batch=="ukb200k" &
                             !(d_ukb_450k.nomiss$eid %in% 
                                 relSample_filter.200k)) 
# UKB250k only
d_ukb250k.nomiss <- subset(d_ukb_450k.nomiss, 
                           d_ukb_450k.nomiss$Batch=="ukb250k" &
                             !(d_ukb_450k.nomiss$eid %in% 
                                 relSample_filter.250k)) 

# related samples 32059: 453375 - 32059 = 421316
d_ukb_450k.nomiss_noRel <- subset(d_ukb_450k.nomiss, 
                            !(d_ukb_450k.nomiss$eid %in% 
                                sam2remove$IID)) 

# total samples: 421,316 for GWAS
#   # N=421316
# fwrite(d_ukb_450k.nomiss_noRel, "/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/CH_phenoCovar.ukb200k_ukb250k_N421316.8Mar2023.tsv.gz",
#        row.names = F, col.names = T, sep="\t", na = 'NA', quote = F)
#   # N=193342
# fwrite(d_ukb200k.nomiss, "/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/CH_phenoCovar.ukb200k_N193342.8Mar2023.tsv.gz",
#        row.names = F, col.names = T, sep="\t", na = 'NA', quote = F)
#   # N=243350
# fwrite(d_ukb250k.nomiss, "/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/CH_phenoCovar.ukb250k_N243350.8Mar2023.tsv.gz",
#        row.names = F, col.names = T, sep="\t", na = 'NA', quote = F)

  ## selected columns to avoid issues in regenie
# fwrite(d_ukb_450k.nomiss_noRel[,c(1,2,31:41,15,16,17,19:29)], "/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/CH_phenoCovar.ukb200k_ukb250k_N421316.27cols.8Mar2023.tsv.gz", 
#        row.names = F, col.names = T, sep="\t", na = 'NA', quote = F)
# fwrite(d_ukb200k.nomiss[,c(1,2,31:41,15,16,17,19:29)], "/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/CH_phenoCovar.ukb200k_N193342.27cols.8Mar2023.tsv.gz",
#        row.names = F, col.names = T, sep="\t", na = 'NA', quote = F)
# fwrite(d_ukb250k.nomiss[,c(1,2,31:41,15,16,17,19:29)], "/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/CH_phenoCovar.ukb250k_N243350.27cols.8Mar2023.tsv.gz",
#        row.names = F, col.names = T, sep="\t", na = 'NA', quote = F)

##################################################

#### VAF GWAS data
ch_var_ukb450.filtered <- subset(ch_var_ukb450, 
                                 (ch_var_ukb450$eid_7089 %in% 
                                    d_ukb_450k.nomiss$FID) & 
                                   !(ch_var_ukb450$eid_7089 %in% 
                                       relSample_filter.VAF_samples) )

ch_var_ukb450.filtered <- merge(d_ukb_450k.nomiss, 
                                ch_var_ukb450.filtered,
                                by.x="eid", by.y="eid_7089" )

# Sort by clone size
ch_var_ukb450.filtered <- ch_var_ukb450.filtered[
  order(ch_var_ukb450.filtered$eid,
        ch_var_ukb450.filtered$AF, 
        decreasing = TRUE), ]

# keep larger clone in individuals with >1 CHIP
ch_var_ukb450.filtered <- ch_var_ukb450.filtered[!duplicated(ch_var_ukb450.filtered$eid),]


sort(table(ch_var_ukb450.filtered$Gene))

ch_var_ukb450.filtered$Gene_Mut <- ifelse(ch_var_ukb450.filtered$Gene=="DNMT3A",
                                     "DNMT3A",
                                     ifelse(ch_var_ukb450.filtered$Gene=="TET2",
                                            "TET2", 
                                            ifelse(ch_var_ukb450.filtered$Gene=="ASXL1",
                                                   "ASXL1",
                                                   ifelse(ch_var_ukb450.filtered$Gene %in% c("PPM1D","TP53"),
                                                          "DDR",
                                                          ifelse(ch_var_ukb450.filtered$Gene %in% c("SF3B1","SRSF2","U2AF1","ZRSR2"),
                                                                 "SF","Other")))))
sort(table(ch_var_ukb450.filtered$Gene_Mut, exclude = NULL), decreasing = T)
# DNMT3A   TET2  Other  ASXL1    DDR     SF 
# 15980   5175   3660   2531   1491    590
names(ch_var_ukb450.filtered)

ch_var_ukb450.filtered <- ch_var_ukb450.filtered[,c(2,3,1,4:29,40:52,54)]
ch_var_ukb450.filtered$varID <- paste(ch_var_ukb450.filtered$CHR, 
                                      ch_var_ukb450.filtered$POS_hg38, 
                                      ch_var_ukb450.filtered$REF, 
                                      ch_var_ukb450.filtered$ALT, 
                                      sep="_")


ukb450k_dnmt3a_qt_vaf <- ch_var_ukb450.filtered[ch_var_ukb450.filtered$Gene=="DNMT3A",]

ukb450k_tet2_qt_vaf   <- ch_var_ukb450.filtered[ch_var_ukb450.filtered$Gene=="TET2",]  

ukb450k_asxl1_qt_vaf   <- ch_var_ukb450.filtered[ch_var_ukb450.filtered$Gene=="ASXL1",]  

ukb450k_dta_qt_vaf   <- ch_var_ukb450.filtered[ch_var_ukb450.filtered$Gene %in% 
                                                 c("DNMT3A" ,"TET2","ASXL1"),]
ukb450k_sf_qt_vaf   <- ch_var_ukb450.filtered[ch_var_ukb450.filtered$Gene %in% 
                                                c("SF3B1","SRSF2","U2AF1","ZRSR2"),]
ukb450k_ddr_qt_vaf   <- ch_var_ukb450.filtered[ch_var_ukb450.filtered$Gene %in% 
                                                c("PPM1D","TP53"),]

# 
# fwrite(ch_var_ukb450.filtered,"/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/ukb450k.CH_qt_vaf.7Mar2023.tsv.gz",
#        row.names=F, col.names=T, sep="\t", quote=F, na="NA")
# fwrite(ukb450k_dnmt3a_qt_vaf,"/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/ukb450k.DNMT3A_qt_vaf.7Mar2023.tsv.gz",
#        row.names=F, col.names=T, sep="\t", quote=F, na="NA")
# fwrite(ukb450k_tet2_qt_vaf,"/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/ukb450k.TET2_qt_vaf.7Mar2023.tsv.gz",
#        row.names=F, col.names=T, sep="\t", quote=F, na="NA")
# fwrite(ukb450k_dta_qt_vaf,"/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/ukb450k.DTA_qt_vaf.7Mar2023.tsv.gz",
#        row.names=F, col.names=T, sep="\t", quote=F, na="NA")
# fwrite(ukb450k_asxl1_qt_vaf,"/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/ukb450k.ASXL1_qt_vaf.7Mar2023.tsv.gz",
#        row.names=F, col.names=T, sep="\t", quote=F, na="NA")
# fwrite(ukb450k_ddr_qt_vaf,"/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/ukb450k.DDR_qt_vaf.7Mar2023.tsv.gz",
#        row.names=F, col.names=T, sep="\t", quote=F, na="NA")
# fwrite(ukb450k_sf_qt_vaf,"/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/ukb450k.SF_qt_vaf.7Mar2023.tsv.gz",
#        row.names=F, col.names=T, sep="\t", quote=F, na="NA")

### All files so far saved in the image
# save.image(file ="/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/CH_phenoCovar.ukb200k_ukb250k.8Mar2023.rda")
load("/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/CH_phenoCovar.ukb200k_ukb250k.8Mar2023.rda")

####################
# WB vs others

#   ## CH 0/1: binary
# ## selected columns to avoid issues in regenie
#   # WB=352058
# fwrite(d_ukb_450k.nomiss_noRel[d_ukb_450k.nomiss_noRel$Ethnic_Background=="WB",c(1,2,31:41,15,16,17,19:29)], 
#        "/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/CH_phenoCovar.WB.ukb200k_ukb250k_N352058.27cols.8Mar2023.tsv.gz",
#        row.names = F, col.names = T, sep="\t", na = 'NA', quote = F)
# # Others=69258
# fwrite(d_ukb_450k.nomiss_noRel[d_ukb_450k.nomiss_noRel$Ethnic_Background=="Others",c(1,2,31:41,15,16,17,19:29)], 
#        "/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/CH_phenoCovar.notWB.ukb200k_ukb250k_N69258.27cols.8Mar2023.tsv.gz",
#        row.names = F, col.names = T, sep="\t", na = 'NA', quote = F)
# 
#   # WB=161151
# fwrite(d_ukb200k.nomiss[d_ukb200k.nomiss$Ethnic_Background=="WB",c(1,2,31:41,15,16,17,19:29)], 
#        "/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/CH_phenoCovar.WB.ukb200k_N161151.27cols.8Mar2023.tsv.gz",
#        row.names = F, col.names = T, sep="\t", na = 'NA', quote = F)
#   # WB=32191
# fwrite(d_ukb200k.nomiss[d_ukb200k.nomiss$Ethnic_Background=="Others",c(1,2,31:41,15,16,17,19:29)], 
#        "/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/CH_phenoCovar.notWB.ukb200k_N32191.27cols.8Mar2023.tsv.gz",
#        row.names = F, col.names = T, sep="\t", na = 'NA', quote = F)
# 
#   # WB=204443
# fwrite(d_ukb250k.nomiss[d_ukb250k.nomiss$Ethnic_Background=="WB",c(1,2,31:41,15,16,17,19:29)], 
#        "/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/CH_phenoCovar.WB.ukb250k_N204443.27cols.8Mar2023.tsv.gz",
#        row.names = F, col.names = T, sep="\t", na = 'NA', quote = F)
# # WB=38907
# fwrite(d_ukb250k.nomiss[d_ukb250k.nomiss$Ethnic_Background=="Others",c(1,2,31:41,15,16,17,19:29)], 
#        "/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/CH_phenoCovar.notWB.ukb250k_N38907.27cols.8Mar2023.tsv.gz",
#        row.names = F, col.names = T, sep="\t", na = 'NA', quote = F)
# 
# ######
# 
#   ## VAF qt
# fwrite(ch_var_ukb450.filtered[ch_var_ukb450.filtered$Ethnic_Background=="WB", 
#                               c(1:3,15,16,31,19:30,40,43,39,44)],
#        "/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/ukb450k.CH_qt_vaf.WB.7Mar2023.tsv.gz",
#        row.names=F, col.names=T, sep="\t", quote=F, na="NA")
# 
# fwrite(ukb450k_dnmt3a_qt_vaf[ukb450k_dnmt3a_qt_vaf$Ethnic_Background=="WB", 
#                              c(1:3,15,16,31,19:30,40,43,39,44)],
#        "/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/ukb450k.DNMT3A_qt_vaf.WB.7Mar2023.tsv.gz",
#        row.names=F, col.names=T, sep="\t", quote=F, na="NA")
# 
# fwrite(ukb450k_tet2_qt_vaf[ukb450k_tet2_qt_vaf$Ethnic_Background=="WB", 
#                            c(1:3,15,16,31,19:30,40,43,39,44)],
#        "/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/ukb450k.TET2_qt_vaf.WB.7Mar2023.tsv.gz",
#        row.names=F, col.names=T, sep="\t", quote=F, na="NA")
# 
# fwrite(ukb450k_dta_qt_vaf[ukb450k_dta_qt_vaf$Ethnic_Background=="WB", 
#                           c(1:3,15,16,31,19:30,40,43,39,44)],
#        "/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/ukb450k.DTA_qt_vaf.WB.7Mar2023.tsv.gz",
#        row.names=F, col.names=T, sep="\t", quote=F, na="NA")
# 
# fwrite(ukb450k_asxl1_qt_vaf[ukb450k_asxl1_qt_vaf$Ethnic_Background=="WB", 
#                             c(1:3,15,16,31,19:30,40,43,39,44)],
#        "/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/ukb450k.ASXL1_qt_vaf.WB.7Mar2023.tsv.gz",
#        row.names=F, col.names=T, sep="\t", quote=F, na="NA")
# 
# fwrite(ukb450k_ddr_qt_vaf[ukb450k_ddr_qt_vaf$Ethnic_Background=="WB", 
#                           c(1:3,15,16,31,19:30,40,43,39,44)],
#        "/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/ukb450k.DDR_qt_vaf.WB.7Mar2023.tsv.gz",
#        row.names=F, col.names=T, sep="\t", quote=F, na="NA")
# 
# fwrite(ukb450k_sf_qt_vaf[ukb450k_sf_qt_vaf$Ethnic_Background=="WB", 
#                          c(1:3,15,16,31,19:30,40,43,39,44)],
#        "/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb450_2023/ukb450k.SF_qt_vaf.WB.7Mar2023.tsv.gz",
#        row.names=F, col.names=T, sep="\t", quote=F, na="NA")



