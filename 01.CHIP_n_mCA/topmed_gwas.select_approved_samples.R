#############
### Select Approved samples
## Topmed
library(data.table)

library("glue")

not_approved <- c("Africa6K:GRU-IRB-PUB-COL-NPU", 
                  "ChildrensHS_GAP:GRU", 
                  "ChildrensHS_IGERA:GRU", 
                  "ChildrensHS_MetaAir:GRU",
                  "GALAI:DS-LD-IRB-COL",
                  "GALAII:DS-LD-IRB-COL",
                  "GENAF:HMB-NPU",
                  "IPF:DS-ILD-IRB-NPU",
                  "IPF:DS-LD-IRB-NPU", 
                  "IPF:DS-PFIB-IRB-NPU", 
                  "IPF:DS-PUL-ILD-IRB-NPU", 
                  "SAGE:DS-LD-IRB-COL",
                  "PCGC_CHD:DS-CHD")

approved_topmed_cohorts <- fread("TOPMed_proposal/parentstudies.csv", 
                                 header=T, sep=",")

approved_status_topmed <- fread("TOPMed_proposal/requested_data.csv", 
                                header=T, sep=",")
## WGS
## with consent annotations
freeze9b_sample_annot <- fread("freeze9b_sample_annot_2021-09-28.txt", header=T)

freeze10a_sample_annot <- fread("freeze.10.draft.sample.mapping.2021dec20.tab.gz", header=T)

sra_topmed <- fread("SRA_ID_mapping/TOPMed_SraRunTable_20210401.txt.gz", header=T)
# PCs 
load("freeze.10b.annotation/AncestryPCs/TOPMed_F10b_PCs.Rdata")
heatmap(cor(PCs[,c(2:21)]))
## 
# load("/Volumes/medpop_esp2/mesbah/datasets/topmed/whi/91227/asp-topmed-dcc/exchange/phs001237_TOPMed_WGS_WHI/Combined_Study_Data/Genotypes/freeze.9b/relatedness/freeze9_pcair_results.RData")


##### select topmed samples
ch_topmed <- fread("TOPMed_CHIP_calls_8-31-2020/TOPMed_CHIPcalls_with_covariates_2020_08_31.tsv", header = T)

varch_topmed <- fread("TOPMed_CHIP_calls_8-31-2020/TOPMed_variant_level_CHIPcalls_with_covariates_2020_08_31.tsv", header = T)

# Josh's filtered phenotype
topmed_ch_anc <- fread("MetaGWAS_N900k/chip_input_01_filter_04_19_2024.tsv", 
                       header = T)

topmed_filtered_ch <- merge(topmed_ch_anc[, c(1:6,17)], 
                            ch_topmed, 
                            by.x="NWD_ID", 
                            by.y="Sample")
summary(topmed_filtered_ch)

## 
fam_topmed <- fread("geno/plink/topmed_N74k.MultiANC.freeze.10b.chr1_22.plink_maf01hwe1e50geno10mind10.fam", header=F)

## 
c("IPF", "PCGC_CHD:DS-CHD")
## c("ChildrensHS_GAP", "ChildrensHS_IGERA", "ChildrensHS_MetaAir")

table(freeze10a_sample_annot$STUDY[freeze10a_sample_annot$CONSENT %in% "DS-LD-IRB-COL"], exclude = NULL)

## 
# samples_apprived <- freeze10a_sample_annot$SAMPLE_ID[paste(freeze10a_sample_annot$STUDY, freeze10a_sample_annot$CONSENT, sep=":") %in% approved_status_topmed$`Parent Study: consent group`[grepl(pattern = "Approved", x = approved_status_topmed$`Review Status`)] ]

samples2exclude <- c(freeze10a_sample_annot$SAMPLE_ID[freeze10a_sample_annot$STUDY  %in% c("ChildrensHS_GAP", "ChildrensHS_IGERA", "ChildrensHS_MetaAir") & freeze10a_sample_annot$CONSENT  %in% c("GRU")],
                     freeze10a_sample_annot$SAMPLE_ID[freeze10a_sample_annot$STUDY  %in% c("GALAI", "GALAII") & freeze10a_sample_annot$CONSENT  %in% c("DS-LD-IRB-COL")],
                     freeze10a_sample_annot$SAMPLE_ID[freeze10a_sample_annot$STUDY  %in% c("GENAF") & freeze10a_sample_annot$CONSENT  %in% c("HMB-NPU")],
                     freeze10a_sample_annot$SAMPLE_ID[freeze10a_sample_annot$STUDY  %in% c("IPF") & freeze10a_sample_annot$CONSENT  %in% c("DS-ILD-IRB-NPU", "DS-LD-IRB-NPU", "DS-PFIB-IRB-NPU", "DS-PUL-ILD-IRB-NPU")], 
                     freeze10a_sample_annot$SAMPLE_ID[freeze10a_sample_annot$STUDY  %in% c("SAGE") & freeze10a_sample_annot$CONSENT  %in% c("DS-LD-IRB-COL")],
                     freeze10a_sample_annot$SAMPLE_ID[freeze10a_sample_annot$STUDY  %in% c("PCGC") & freeze10a_sample_annot$CONSENT  %in% c("DS-CHD")]
)
## exclude 
table(topmed_filtered_ch$NWD_ID %in% samples2exclude)
# FALSE  TRUE 
# 74622   344

### final Topmed data frame
# exclude 344 samples w/o approval
topmed_filtered_ch <- merge(topmed_ch_anc[!(topmed_ch_anc$NWD_ID %in% samples2exclude), 
                                          c(1:6,17)], 
                            freeze10a_sample_annot, 
                            by.x="NWD_ID", 
                            by.y="SAMPLE_ID")

table(topmed_filtered_ch$INFERRED_SEX, topmed_filtered_ch$SEX, 
      exclude = NULL)

table(freeze9b_sample_annot$sex[freeze9b_sample_annot$sample.id %in% topmed_filtered_ch$NWD_ID[topmed_filtered_ch$INFERRED_SEX==1]])

table(freeze9b_sample_annot$sex[freeze9b_sample_annot$sample.id %in% topmed_filtered_ch$NWD_ID[topmed_filtered_ch$INFERRED_SEX==2]])

## exclude Sex mismatch 
table(topmed_filtered_ch$INFERRED_SEX==1 & topmed_filtered_ch$SEX=="male")

topmed_filtered_ch$Genetic_Sex <- ifelse(topmed_filtered_ch$INFERRED_SEX==1 & 
                                           topmed_filtered_ch$SEX=="male", 
                                         "Male", 
                                         ifelse(topmed_filtered_ch$INFERRED_SEX==2 & 
                                                  topmed_filtered_ch$SEX=="female", 
                                                "Female", 
                                                "NA"))

table(topmed_filtered_ch$Genetic_Sex, exclude = NULL)
# Female   Male     NA 
# 43885  29318   1419
N=73203
Female=43885 # 
Male=29318 # 40.05027%

table(topmed_filtered_ch$knn[!is.na(topmed_filtered_ch$Genetic_Sex)], exclude = NULL)
# AFR   AMR   EAS   EUR   SAS 
# 20371  8087  4664 41146   354

round(prop.table(table(topmed_filtered_ch$knn[!is.na(topmed_filtered_ch$Genetic_Sex)], exclude = NULL)) *100, 2)

# AFR   AMR   EAS   EUR   SAS 
# 27.30 10.84  6.25 55.14  0.47 

### Annotate CHIP

sort(table(varch_topmed$Gene[varch_topmed$Sample %in% topmed_filtered_ch$NWD_ID], exclude = NULL))

# min VAF in TOPMed >=5%
topmed_filtered_ch$hasCHvaf02 <- ifelse(topmed_filtered_ch$NWD_ID %in% varch_topmed$Sample, 1,0) 

table(topmed_filtered_ch$hasCHvaf02, topmed_filtered_ch$hasCHIP, 
      exclude=NULL)

# VAF>=10%
topmed_filtered_ch$hasCHvaf10 <- ifelse(topmed_filtered_ch$NWD_ID %in% varch_topmed$Sample[varch_topmed$VAF>=0.1], 
                                        1, 
                                        ifelse(topmed_filtered_ch$hasCHvaf02==0, 
                                               0 , NA)) 

# 0     1
# 0 69951     0
# 1     0  4671
table(topmed_filtered_ch$hasCHvaf10, exclude=NULL)
# 0     1  <NA> 
#  69951  4074   597

# DNMT3A
topmed_filtered_ch$hasDNMT3A <- ifelse(topmed_filtered_ch$NWD_ID %in% varch_topmed$Sample[varch_topmed$Gene=="DNMT3A"], 
                                       1, 
                                       ifelse(topmed_filtered_ch$hasCHvaf02==0, 
                                              0 , NA)) 

table(topmed_filtered_ch$hasDNMT3A, exclude=NULL)
# 0     1  <NA> 
#  69951  2448  2223

# TET2
topmed_filtered_ch$hasTET2 <- ifelse(topmed_filtered_ch$NWD_ID %in% varch_topmed$Sample[varch_topmed$Gene == "TET2"], 
                                     1, 
                                     ifelse(topmed_filtered_ch$hasCHvaf02==0, 
                                            0 , NA)) 


table(topmed_filtered_ch$hasTET2, exclude=NULL)
# 0     1  <NA> 
#  69951   890  3781
# ASXL1
topmed_filtered_ch$hasASXL1 <- ifelse(topmed_filtered_ch$NWD_ID %in% varch_topmed$Sample[varch_topmed$Gene == "ASXL1"], 
                                      1, 
                                      ifelse(topmed_filtered_ch$hasCHvaf02==0, 
                                             0 , NA)) 


table(topmed_filtered_ch$hasASXL1, exclude=NULL)
# 0     1  <NA> 
#  69951   363  4308

# Splicing Factors
topmed_filtered_ch$hasSF <- ifelse(topmed_filtered_ch$NWD_ID %in% 
                                     varch_topmed$Sample[varch_topmed$Gene %in% 
                                                           c("SF3B1", "SRSF2", "U2AF1", "ZRSR2")], 
                                   1, 
                                   ifelse(topmed_filtered_ch$hasCHvaf02==0, 
                                          0 , NA)) 


table(topmed_filtered_ch$hasSF, exclude=NULL)
# 0     1  <NA> 
#  69951   266  4405

## DNA Damage Genes: c("PPM1D", "TP53")
topmed_filtered_ch$hasDDR <- ifelse(topmed_filtered_ch$NWD_ID %in% varch_topmed$Sample[varch_topmed$Gene %in% c("PPM1D", "TP53")], 
                                    1, 
                                    ifelse(topmed_filtered_ch$hasCHvaf02==0, 
                                           0 , NA))  

table(topmed_filtered_ch$hasDDR, exclude=NULL)
##     0     1  <NA> 
# 69951   262  4409 

table(topmed_filtered_ch$knn[!is.na(topmed_filtered_ch$Genetic_Sex)], 
      topmed_filtered_ch$hasDDR[!is.na(topmed_filtered_ch$Genetic_Sex)])

topmed_filtered_ch %>% filter(!is.na(Genetic_Sex)) %>% group_by(knn) %>% summarise(N= n())
# A tibble: 5 × 2
# knn       N
# <chr> <int>
#   1 AFR   20371
# 2 AMR    8087
# 3 EAS    4664
# 4 EUR   41146
# 5 SAS     354
topmed_filtered_ch %>% filter(!is.na(Genetic_Sex)) %>% group_by(knn) %>% tally(hasDDR)


## Age at blood draw
topmed_filtered_ch <- merge(topmed_filtered_ch, chip[, c(1,6)], 
                            by="NWD_ID")

topmed_filtered_ch$AGE=NULL

topmed_filtered_ch$AgeAtBloodDraw <- round(topmed_filtered_ch$AgeAtBloodDraw, 2)

topmed_filtered_ch$Sqrd_AgeAtBloodDraw <- round(topmed_filtered_ch$AgeAtBloodDraw^2,2) 
# TOPMED PCs
PCs$NWD_ID <- as.character(PCs$SAMPLE_ID)

topmed_filtered_ch <- merge(topmed_filtered_ch, 
                            PCs[,c(22, 2:21)], 
                            by = "NWD_ID")

## Save file
date_label = format(Sys.Date(), "%m_%d_%Y")

data.table::fwrite(topmed_filtered_ch, 
                   glue("pheno/topmed_input_gwas_{date}.tsv.gz", date = date_label), 
                   row.names = F, 
                   col.names = T, 
                   sep="\t", 
                   quote = F, 
                   na = "NA")



# topmed_filtered_ch <- data.table::fread(glue("pheno/topmed_input_gwas_{date}.tsv.gz", date = date_label))
##############################
topmed_filtered_ch <- data.table::fread("pheno/topmed_input_gwas_04_23_2024.tsv.gz", header=T)
summary(topmed_filtered_ch$AgeAtBloodDraw)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 36.00   51.00   60.00   59.99   69.00   98.00 
round(sd(topmed_filtered_ch$AgeAtBloodDraw),4)
# 11.1895
table(topmed_filtered_ch$Genetic_Sex, exclude = NULL)
# Female   Male   <NA> 
#   43885  29318   1419
round(prop.table(table(topmed_filtered_ch$Genetic_Sex, exclude = NULL))*100,1)
# Female   Male   <NA> 
#   58.8   39.3    1.9
table(topmed_filtered_ch$knn, exclude = NULL)
# AFR   AMR   EAS   EUR   SAS 
# 20371  8087  4664 41146   354
topmed_filtered_ch.v2 <- data.table::fread("pheno/topmed_input_gwas_04_23_2024.v02.tsv.gz", header=T)
summary(topmed_filtered_ch.v2$AgeAtBloodDraw)

## CHIP Genes in TOPMed GWAS samples
table(unique(varch_topmed$Sample)%in% topmed_filtered_ch$NWD_ID)
# FALSE  TRUE 
# 831  4671 
nrow(varch_topmed)
qcd.topmed_var_table <- subset(varch_topmed, 
                               varch_topmed$Sample%in% topmed_filtered_ch$NWD_ID)
nrow(qcd.topmed_var_table)

table(table(qcd.topmed_var_table$Sample))
head(sort(table(qcd.topmed_var_table$Gene), decreasing = T),20)

head(sort(round(prop.table(table(qcd.topmed_var_table$Gene))*100,2), decreasing = T),20)

##############################
## TOPMed data in GWAS #####
topmed_supl_table1 <- topmed_filtered_ch %>%
  group_by(STUDY.y) %>%
  summarise(
    N       = n(),
    avgAge  = round(mean(AgeAtBloodDraw, na.rm = TRUE),1),
    sd_Age  = round(sd(AgeAtBloodDraw, na.rm = TRUE),1),
    Seq_Ctr   = paste(unique(SEQ_CTR),   collapse = ", "),
    subStud  = paste(unique(STUDY.x),   collapse = ", "),
    PHS = paste(unique(PHS),   collapse = ", "),
    CONSENT = paste(unique(CONSENT),   collapse = ", "),
    AFR = sum(knn=="AFR"),
    AMR = sum(knn=="AMR"),
    EAS = sum(knn=="EAS"),
    EUR = sum(knn=="EUR"),
    SAS = sum(knn=="SAS")
  )

# write.csv(topmed_supl_table1, "../../../NG_v2/Figs/chip_gwas_topmedData_supl_table1.csv", row.names = F)

############################
##
topmed_filtered_ch$FID <- 0

topmed_filtered_ch$IID <- topmed_filtered_ch$NWD_ID


names(topmed_filtered_ch)

### Sequencing Centers
table(topmed_filtered_ch$SEQ_CTR, exclude = NULL)
# Baylor    Broad Illumina Macrogen     NYGC       UW    WashU 
# 23984    31370      714      854     1726    11360     4614 
topmed_filtered_ch$Sequencing_Center <- ifelse(!(topmed_filtered_ch$SEQ_CTR %in% c("Baylor", "Broad","UW") ), "Others", topmed_filtered_ch$SEQ_CTR)

table(topmed_filtered_ch$Sequencing_Center, exclude = NULL)

## TOPMed Phase
table(topmed_filtered_ch$PHASE, exclude = NULL)

topmed_filtered_ch$TOPMed_Phase <- ifelse(!(topmed_filtered_ch$PHASE %in% c(1,2,3) ), "Others", topmed_filtered_ch$PHASE)

table(topmed_filtered_ch$TOPMed_Phase, exclude = NULL)

### TOPMed Study
# N>3000
n3k_studies <- c(names(table(topmed_filtered_ch$STUDY.y)[(table(topmed_filtered_ch$STUDY.y)>=3000)]))

topmed_filtered_ch$STUDY <- ifelse( !(topmed_filtered_ch$STUDY.y %in% n3k_studies), "not_3k",  topmed_filtered_ch$STUDY.y)

table(topmed_filtered_ch$STUDY, exclude = NULL)

## MultiAnc
# fwrite(topmed_filtered_ch[, c(55, 56, 1, 26:34, 25, 6, 57:59, 35:54)], 
#                    glue("pheno/TOPMed_GWAS_pheno.MultiANC.{date}.tsv.gz", date = date_label), 
#                    row.names = F, 
#                    col.names = T, 
#                    sep="\t", 
#                    quote = F, 
#                    na = "NA")
### MultiANC

topmed_filtered_ch$GenANC <- ifelse(topmed_filtered_ch$knn %in% c("SAS", "EAS"), "Others", topmed_filtered_ch$knn) 

table(topmed_filtered_ch$GenANC, exclude =NULL)
## 
vars <- "v02"

ANC="MultiANC"
fwrite(topmed_filtered_ch[!is.na(topmed_filtered_ch$Genetic_Sex), 
                          c(55, 56, 1, 26:34, 25, 60, 57:59, 35:54)], 
       glue("pheno/TOPMed_GWAS_pheno.{ANC}.{date}.{vars}.tsv.gz", date = date_label), 
       row.names = F, 
       col.names = T, 
       sep="\t", 
       quote = F, 
       na = "NA")

## Female
ANC="Female"
fwrite(topmed_filtered_ch[topmed_filtered_ch$Genetic_Sex==ANC, 
                          c(55, 56, 1, 26:34, 25, 60, 57:59, 35:54)], 
       glue("pheno/TOPMed_GWAS_pheno.{ANC}.{date}.{vars}.tsv.gz", date = date_label), 
       row.names = F, 
       col.names = T, 
       sep="\t", 
       quote = F, 
       na = "NA")

## Male
ANC="Male"
fwrite(topmed_filtered_ch[topmed_filtered_ch$Genetic_Sex==ANC, c(55, 56, 1, 26:34, 25, 60, 57:59, 35:54)], 
       glue("pheno/TOPMed_GWAS_pheno.{ANC}.{date}.{vars}.tsv.gz", date = date_label), 
       row.names = F, 
       col.names = T, 
       sep="\t", 
       quote = F, 
       na = "NA")


## AFR
ANC="AFR"
fwrite(topmed_filtered_ch[!is.na(topmed_filtered_ch$Genetic_Sex) & (topmed_filtered_ch$knn==ANC), 
                          c(55, 56, 1, 26:34, 25, 60, 57:59, 35:54)], 
       glue("pheno/TOPMed_GWAS_pheno.{ANC}.{date}.{vars}.tsv.gz", date = date_label), 
       row.names = F, 
       col.names = T, 
       sep="\t", 
       quote = F, 
       na = "NA")

## AMR
ANC="AMR"
fwrite(topmed_filtered_ch[!is.na(topmed_filtered_ch$Genetic_Sex) & (topmed_filtered_ch$knn==ANC), 
                          c(55, 56, 1, 26:34, 25, 60, 57:59, 35:54)], 
       glue("pheno/TOPMed_GWAS_pheno.{ANC}.{date}.{vars}.tsv.gz", date = date_label), 
       row.names = F, 
       col.names = T, 
       sep="\t", 
       quote = F, 
       na = "NA")

## EUR
ANC="EUR"
fwrite(topmed_filtered_ch[!is.na(topmed_filtered_ch$Genetic_Sex) & (topmed_filtered_ch$knn==ANC), 
                          c(55, 56, 1, 26:34, 25, 60, 57:59, 35:54)], 
       glue("pheno/TOPMed_GWAS_pheno.{ANC}.{date}.{vars}.tsv.gz", date = date_label), 
       row.names = F, 
       col.names = T, 
       sep="\t", 
       quote = F, 
       na = "NA")

# ## SAS
# ANC="SAS"
# fwrite(topmed_filtered_ch[!is.na(topmed_filtered_ch$Genetic_Sex) & (topmed_filtered_ch$knn==ANC), 
#                           c(55, 56, 1, 26:34, 25, 60, 57:59, 35:54)], 
#        glue("pheno/TOPMed_GWAS_pheno.{ANC}.{date}.{vars}.tsv.gz", date = date_label), 
#        row.names = F, 
#        col.names = T, 
#        sep="\t", 
#        quote = F, 
#        na = "NA")
# 
# ## EAS
# ANC="EAS"
# fwrite(topmed_filtered_ch[!is.na(topmed_filtered_ch$Genetic_Sex) & (topmed_filtered_ch$knn==ANC), 
#                           c(55, 56, 1, 26:34, 25, 60, 57:59, 35:54)], 
#        glue("pheno/TOPMed_GWAS_pheno.{ANC}.{date}.{vars}.tsv.gz", date = date_label), 
#        row.names = F, 
#        col.names = T, 
#        sep="\t", 
#        quote = F, 
#        na = "NA")
# 
##########################################


##############

# AFR   AMR   EAS   EUR   SAS 
# 20385  8111  4664 41452   354 
round(prop.table(table(topmed_ch_anc$knn))*100,2)
# AFR   AMR   EAS   EUR   SAS 
# 27.19 10.82  6.22 55.29  0.4
