# use R-4.1
library(data.table)

	# p21022 = Age at recruitment, 
	# p31= Sex, 
	# p22001 = Genetic sex
	# p22019 = Sex chromosome aneuploidy
	# p22006 = Genetic ethnic grouping "Caucasian"
	# p21000_i0 = Ethnic background 
	# p21001_i0 = BMI
	# p23104_i0 = Body mass index (BMI)
	# p20160_i0 = Ever smoked
	# p20116_i0=Smoking status 
	# p22000 = Genotype measurement batch BiLEVE vs Axiom
	# from RAP
d_base <- fread("/medpop/esp2/mesbah/projects/gxe/UKB_baseline/UKB500k_baseline_info_recoded_participant.tsv.gz")

	# Genotype array Batch
d_base$GenoBatch <- ifelse(grepl(pattern="UKBiLEVEAX", x=d_base$p22000),"UKBiLEVEAX","Axiom")
# table(d_base$GenoBatch)
#     Axiom UKBiLEVEAX 
#   452476      49934
	# Age
d_base$Age_at_recruitment <- d_base$p21022
# summary(d_base$Age_at_recruitment)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  37.00   50.00   58.00   56.53   63.00   73.00       1

	# Sex
d_base$Genetic_Sex[d_base$p22001==d_base$p31] <- d_base$p22001[d_base$p22001==d_base$p31]
# exclude sex mismatch
#table(!is.na(d_base$Genetic_Sex))
# FALSE   TRUE 
# 14611 487799

	# Sex chromosome aneuploidy
d_base$Aneuploidy <- ifelse(d_base$p22019=="Yes", "Aneuploidy","No_Aneuploidy")
# table(d_base$Aneuploidy)
#   Aneuploidy No_Aneuploidy 
#          651        501759 

	
#################
### load big csv file
d_pcs <- fread("/medpop/esp2/projects/UK_Biobank/baskets/2008463/ukb47823.csv.gz", select = c('eid', '22009-0.1', '22009-0.2', '22009-0.3', '22009-0.4', '22009-0.5', '22009-0.6', '22009-0.7', '22009-0.8', '22009-0.9', '22009-0.10'))

names(d_pcs) <- c("eid", paste0("PC",1:10))

## 
d_ukb_all <- merge(d_base, d_pcs, by="eid")

# save.image(file="/medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/ukb500k_covariates.rda")
# 
# fwrite(d_ukb_all,"/medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/ukb_pheno_covar.tsv", row.names=F, col.names=T, sep="\t", quote=F, na="NA")

############################### Annotate CHIP
# Load UKB data
# "d_base"    "d_pcs"     "d_ukb_all"
load("/medpop/esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/ukb500k_covariates.rda")
# load("~/Documents/Project/CHIP_GWAS/rerun/ukb450k/ukb500k_covariates.rda")

# VAF>=2%
## ukb450k no NA
ukb450_sampleids <- fread("/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/450k/ukb450k.eid.list", 
                          header = F)
table(d_ukb_all$eid %in% ukb450_sampleids$V1)
# FALSE   TRUE 
# 47618 454792 
table(d_ukb_all$eid %in% ukb450_sampleids$V1 & !is.na(d_ukb_all$PC1))
# FALSE   TRUE 
# 48011 454399
d_ukb_450k <- subset(d_ukb_all, 
                     d_ukb_all$eid %in% ukb450_sampleids$V1)
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
var_chip.ukb450.CHIPonly <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/chipVar_CV_AB.ukb450k_all.csv.gz")

head(sort(table(c(var_chip.ukb450.CHIPonly$Gene.refGene) ), decreasing = T),15)
# DNMT3A   TET2  ASXL1  PPM1D   TP53  SRSF2   JAK2  ASXL2  SF3B1    NF1 
# 16009   6059   2915    949    812    351    288    266    263    216 
# GNB1    CBL CREBBP   GNAS   CUX1 
# 180    179    170    161    143
## CHIP
d_ukb_450k$hasCHIP <- ifelse(d_ukb_450k$eid %in% 
                               unique(var_chip.ukb450.CHIPonly$Broad_ID),1,
                             ifelse(!(d_ukb_450k$eid %in% 
                                        unique(var_chip.ukb450.CHIPonly$Broad_ID)),
                                    0,NA))

(table(d_ukb_450k$hasCHIP, exclude = NULL))
# 0      1 
# 425665  29127 
prop.table(table(d_ukb_450k$hasCHIP, exclude = NULL))
# 0          1 
# 0.93595534 0.06404466

# DNMT3A
d_ukb_450k$hasDNMT3A <- ifelse(d_ukb_450k$eid %in% 
                                 var_chip.ukb450.CHIPonly$Broad_ID[var_chip.ukb450.CHIPonly$Gene.refGene=="DNMT3A"], 1, 
                               ifelse(d_ukb_450k$hasCHIP==0,0,NA))
(table(d_ukb_450k$hasDNMT3A, exclude = NULL))
# # 
# 0      1   <NA> 
#   425665  15555  13572
# TET2
d_ukb_450k$hasTET2 <- ifelse(d_ukb_450k$eid %in% 
                               var_chip.ukb450.CHIPonly$Broad_ID[var_chip.ukb450.CHIPonly$Gene.refGene=="TET2"], 1, 
                             ifelse(d_ukb_450k$hasCHIP==0,0,NA))
(table(d_ukb_450k$hasTET2, exclude = NULL))
#   0     1  <NA> 
#   425665   5871  23256 

# ASXL1
d_ukb_450k$hasASXL1 <- ifelse(d_ukb_450k$eid %in% 
                                var_chip.ukb450.CHIPonly$Broad_ID[var_chip.ukb450.CHIPonly$Gene.refGene=="ASXL1"], 1, 
                              ifelse(d_ukb_450k$hasCHIP==0,0,NA))
(table(d_ukb_450k$hasASXL1, exclude = NULL))
# 0     1  <NA> 
#   425665   2892  26235
# PPM1D
d_ukb_450k$hasPPM1D <- ifelse(d_ukb_450k$eid %in% 
                                var_chip.ukb450.CHIPonly$Broad_ID[var_chip.ukb450.CHIPonly$Gene.refGene=="PPM1D"], 1, 
                              ifelse(d_ukb_450k$hasCHIP==0,0,NA))
(table(d_ukb_450k$hasPPM1D, exclude = NULL))
# 0     1  <NA> 
#   425665    933  28194
# TP53
d_ukb_450k$hasTP53 <- ifelse(d_ukb_450k$eid %in% 
                               var_chip.ukb450.CHIPonly$Broad_ID[var_chip.ukb450.CHIPonly$Gene.refGene=="TP53"], 1, 
                             ifelse(d_ukb_450k$hasCHIP==0,0,NA))
(table(d_ukb_450k$hasTP53, exclude = NULL))
# 0     1  <NA> 
#   425665    808  28319
# SRSF2
d_ukb_450k$hasSRSF2 <- ifelse(d_ukb_450k$eid %in% 
                                var_chip.ukb450.CHIPonly$Broad_ID[var_chip.ukb450.CHIPonly$Gene.refGene=="SRSF2"], 1, 
                              ifelse(d_ukb_450k$hasCHIP==0,0,NA))
(table(d_ukb_450k$hasSRSF2, exclude = NULL))
# 0     1  <NA> 
#   425665    350  28777
# JAK2
d_ukb_450k$hasJAK2 <- ifelse(d_ukb_450k$eid %in% 
                               var_chip.ukb450.CHIPonly$Broad_ID[var_chip.ukb450.CHIPonly$Gene.refGene=="JAK2"], 1, 
                             ifelse(d_ukb_450k$hasCHIP==0,0,NA))
(table(d_ukb_450k$hasJAK2, exclude = NULL))
# 0     1  <NA> 
#   425665    288  28839 
# ASXL2
d_ukb_450k$hasASXL2 <- ifelse(d_ukb_450k$eid %in% 
                                var_chip.ukb450.CHIPonly$Broad_ID[var_chip.ukb450.CHIPonly$Gene.refGene=="ASXL2"], 1, 
                              ifelse(d_ukb_450k$hasCHIP==0,0,NA))
(table(d_ukb_450k$hasASXL2, exclude = NULL))
# 0     1  <NA> 
#   425665    264  28863 
# SF3B1
d_ukb_450k$hasSF3B1 <- ifelse(d_ukb_450k$eid %in% 
                                var_chip.ukb450.CHIPonly$Broad_ID[var_chip.ukb450.CHIPonly$Gene.refGene=="SF3B1"], 1, 
                              ifelse(d_ukb_450k$hasCHIP==0,0,NA))
(table(d_ukb_450k$hasSF3B1, exclude = NULL))
# 0     1  <NA> 
#   425665    262  28865 
# NF1
d_ukb_450k$hasNF1 <- ifelse(d_ukb_450k$eid %in% 
                              var_chip.ukb450.CHIPonly$Broad_ID[var_chip.ukb450.CHIPonly$Gene.refGene=="NF1"], 1, 
                            ifelse(d_ukb_450k$hasCHIP==0,0,NA))
(table(d_ukb_450k$hasNF1, exclude = NULL))
# 0     1  <NA> 
#   425665    212  28915
# GNB1
d_ukb_450k$hasGNB1 <- ifelse(d_ukb_450k$eid %in% 
                               var_chip.ukb450.CHIPonly$Broad_ID[var_chip.ukb450.CHIPonly$Gene.refGene=="GNB1"], 1, 
                             ifelse(d_ukb_450k$hasCHIP==0,0,NA))
(table(d_ukb_450k$hasGNB1, exclude = NULL))
# 0     1  <NA> 
#   425665    180  28947
# CBL
d_ukb_450k$hasCBL <- ifelse(d_ukb_450k$eid %in% 
                              var_chip.ukb450.CHIPonly$Broad_ID[var_chip.ukb450.CHIPonly$Gene.refGene=="CBL"], 1, 
                            ifelse(d_ukb_450k$hasCHIP==0,0,NA))
(table(d_ukb_450k$hasCBL, exclude = NULL))
# 0     1  <NA> 
#   425665    179  28948
### 
d_ukb_450k$FID <- d_ukb_450k$eid
d_ukb_450k$IID <- d_ukb_450k$eid
d_ukb_450k <- d_ukb_450k[, c(40,41,1:39)]

# Ancestry: British-Caucasian coded as "WB" (White British) vs. "Others"  
d_ukb_450k$Ethnic_Background <- ifelse(d_ukb_450k$self_ethnic.p21000_i0=="British" & d_ukb_450k$Genetic_ethnic.p22006=="Caucasian", "WB", "Others")
# Age^2
d_ukb_450k$sqrtAge_at_recruitment <- d_ukb_450k$Age_at_recruitment^2
str(d_ukb_450k)
## 0 sample withdrawn 
ukb_withdrawn <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/w7089_20220222.csv")
table(d_ukb_450k$eid %in% ukb_withdrawn$V1, exclude = NULL)
# FALSE 
# 454792
# fwrite(d_ukb_450k, "/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/chipCall_CV_AB.ukb450k_all.tsv.gz",
#        row.names = F, col.names = T, sep="\t", na = 'NA', quote = F)

#################################################### 
# d_ukb_450k <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/chipCall_CV_AB.ukb450k_all.tsv.gz")
### Step by step sample filtering for GWAS
## UKB450k samples N=454792
# no PCs 393: 454792 - 393 = 454399
d_ukb_450k.nomiss <- subset(d_ukb_450k, 
                            !is.na(d_ukb_450k$PC1) )
# Sex mismatch 94: 454399 - 94 = 454305
d_ukb_450k.nomiss <- subset(d_ukb_450k.nomiss, 
                            !is.na(d_ukb_450k.nomiss$Genetic_Sex))
## Aneuploidy 430: 454305 - 430 = 453875
d_ukb_450k.nomiss <- subset(d_ukb_450k.nomiss, 
                            d_ukb_450k.nomiss$Aneuploidy == "No_Aneuploidy" )
# withdrawn 0: 453875 - 0 = 453875
ukb_withdrawn <- fread("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/w7089_20220222.csv")
table(d_ukb_450k.nomiss$eid %in% ukb_withdrawn$V1)

## Exclude 1st/2nd/3rd degree related samples
geno_related <- fread("~/Documents/Project/CHIP_GWAS/PhenoCov_prep/ukb7089_rel_s488264.dat", 
                      header = T, stringsAsFactors = F)

library(plinkQC)
# related samples 32108: 453875 - 32108 = 421767
geno_related_sam <- relatednessFilter(geno_related, 
                                      relatednessTh = 0.0884,
                                      relatednessIID1 = "ID1",
                                      relatednessIID2 = "ID2",
                                      relatednessRelatedness = "Kinship",
                                      verbose = FALSE)
sam2remove <- geno_related_sam$failIDs
table(d_ukb_450k.nomiss$eid %in% sam2remove$IID)
# FALSE   TRUE 
# 421767  32108
d_ukb_450k.nomiss <- subset(d_ukb_450k.nomiss, !(d_ukb_450k.nomiss$eid %in% sam2remove$IID)) 

d_ukb_450k.nomiss$sqrtAge_at_recruitment <- d_ukb_450k.nomiss$Age_at_recruitment^2
str(d_ukb_450k.nomiss)
# total samples: 421,767 for GWAS
fwrite(d_ukb_450k.nomiss, "/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/chipCallnomiss_CV_AB.ukb450k_all.csv.gz",
       row.names = F, col.names = T, sep=",", na = 'NA', quote = F)
  # Tab-separated file for GWAS: 
# ids1-3,sex,age,age2,ethnicity,batch, pc1-10, chip
d_ukb450k_gwas <- d_ukb_450k.nomiss[,c(1:3,17,16,42,43,15,19:41)]
names(d_ukb450k_gwas) <- c(names(d_ukb450k_gwas)[1:3],
                           "Gender", "Age", "Age2", 
                           "Ethnicity", names(d_ukb450k_gwas)[8:31])
prop.table(table(d_ukb450k_gwas$Gender))*100
# Female     Male 
# 53.93665 46.06335
prop.table(table(d_ukb450k_gwas$Ethnicity))*100
# Others       WB 
# 16.45956 83.54044
table(d_ukb450k_gwas$GenoBatch)
# Axiom UKBiLEVEAX 
# 383344      38423
fwrite(d_ukb450k_gwas, "/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/chipCallnomiss_CV_AB.ukb450k_all.tsv.gz",
       row.names = F, col.names = T, sep="\t", na = 'NA', quote = F)
# save.image(file ="/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/ukb450.22Jul2022.all_files.rda")
# load("/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/ukb450.22Jul2022.all_files.rda")
# ukb200k_eid <- fread("/Volumes/medpop_esp2/mesbah/datasets/CHIP/UKBB/ukb200k_eid.csv.gz")
# d_ukb200k_gwas <- subset(d_ukb450k_gwas, 
#                          (d_ukb450k_gwas$eid %in% ukb200k_eid$eid_7089) )
# fwrite(d_ukb200k_gwas, "/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/chipCallnomiss_CV_AB.ukb200k_all.tsv.gz",
#        row.names = F, col.names = T, sep="\t", na = 'NA', quote = F)
# 
# 
# d_ukb250k_gwas <- subset(d_ukb450k_gwas, 
#                          !(d_ukb450k_gwas$eid %in% ukb200k_eid$eid_7089) )
# fwrite(d_ukb250k_gwas, "/Volumes/medpop_esp2/mesbah/projects/Meta_GWAS/rerun/ukb450k/chipCallnomiss_CV_AB.ukb250k_all.tsv.gz",
#        row.names = F, col.names = T, sep="\t", na = 'NA', quote = F)
# 

###################################################################################################


