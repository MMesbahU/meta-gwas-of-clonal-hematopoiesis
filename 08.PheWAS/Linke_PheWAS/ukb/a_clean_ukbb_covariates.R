###### 01 - library ######
library(dplyr)
library(data.table)

gc()
rm(list=ls())

setwd("/medpop/esp2/lli/chip_gwas_rev/")



##### 02 - files #####
covs <- fread("/medpop/esp2/mzekavat/UKBB/ukbb_PhenoFile.ALL_500k.UpdatedIncdPhenos_202020.REAL.txt")
covs$IDs_toRemove_SampleQCv3 = ifelse((!(covs$Submitted_Gender == covs$Inferred_Gender) |  covs$Non_Consented== 1),1,0)
covs = covs[which(covs$IDs_toRemove_SampleQCv3==0),]

grep("townsend", names(covs), ignore.case = T, value=T)

covs_short <- covs %>%
    select('id', 'smok_detailed_','age', 'age2', 
           'Sex_numeric','SmokingStatusv2','townsend', 
           'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 
           'PC7', 'PC8', 'PC9', 'PC10', "Prev_Maryam_HemeCa_Phenos")

write.table(covs_short, "Data/ukb_phewas/ukb_phewas_covs.tsv", sep="\t", row.names = F, quote = F)



## for cytopenia
blood_counts_fields <- c("Haemoglobin.concentration", "Platelet.count", "Neutrophill.count")
covs_bloodcounts <- covs %>%
    select(id, all_of(blood_counts_fields))

write.table(covs_bloodcounts, "Data/ukb_phewas/ukb_phewas_bloodcounts.tsv", sep="\t", row.names = F, quote = F)



### all of the ukbb phecodes data
all_phecodes <- list.files("/medpop/esp2/mzekavat/UKBB/PhenoFiles/PheCODEs_new/")
phecodes_cat <- all_phecodes[!grepl("2021-01-08_ukb_phecode_(incident|any|prevalent)", all_phecodes)]
# "/medpop/esp2/mzekavat/UKBB/PhenoFiles/PheCODEs_new/2021-01-08_ukb_phecode_incident_March2020fu.csv" 
# is coded binary for incidence cases
# incident.days is for number of days till incidence
# 2021-01-08_ukb_phecode_any_March2020fu.csv also binary

phecodes_cat_final <- gsub("2021-01-08_ukb_phecode_|_March2020fu.csv", "", phecodes_cat)

write.table(phecodes_cat_final, "Data/ukb_phewas/phecode_list.tsv", sep="\t", row.names = F, quote = F, col.names = F)





##### 03 - cleaning covariates data #####
covs <- fread("Data/ukb_phewas/ukb_phewas_covs.tsv", data.table = F)
covs_bloodcounts <- fread("Data/ukb_phewas/ukb_phewas_bloodcounts.tsv", data.table = F)

# no unit conversion is required
covs_bloodcounts1 <- covs_bloodcounts %>%
    left_join(covs %>% select(id, Sex_numeric), by = "id") %>%
    mutate(
        anemia = case_when(
            # Anemia: Hgb <12 g/dL for females (Sex_numeric == 0), <13 g/dL for males (Sex_numeric == 1)
            is.na(Haemoglobin.concentration) ~ NA,
            (Sex_numeric == 0 & Haemoglobin.concentration < 12) | (Sex_numeric == 1 & Haemoglobin.concentration < 13) ~ 1,
            (Sex_numeric == 0 & Haemoglobin.concentration >= 12) | (Sex_numeric == 1 & Haemoglobin.concentration >= 13) ~ 0,
            # when sex is NA
            Haemoglobin.concentration >= 13 ~ 0,
            Haemoglobin.concentration < 12  ~ 1,
            TRUE ~ NA
        ),
        # Thrombocytopenia: platelet count < 150 K/uL
        thrombocytopenia = case_when(
            Platelet.count < 150 ~ 1, 
            Platelet.count >= 150 ~ 0, 
            TRUE ~ NA
        ),
        # Neutropenia: ANC < 1.8 K/uL
        neutropenia = case_when(
            Neutrophill.count < 1.8 ~ 1,
            Neutrophill.count >= 1.8 ~ 0,
            TRUE ~ NA
        ),
        cytopenia = case_when(
            anemia == 0 & thrombocytopenia == 0 & neutropenia == 0 ~ 0,
            anemia == 1 | thrombocytopenia == 1 | neutropenia == 1 ~ 1,
            TRUE ~ NA
        )
    )
    
table(covs_bloodcounts1$anemia, useNA = "ifany") # 14983 NAs
summary(covs_bloodcounts1$Haemoglobin.concentration)  # 14983 NAs

table(covs_bloodcounts1$thrombocytopenia, useNA = "ifany")
summary(covs_bloodcounts1$Platelet.count)  # 14986 NAs


table(covs_bloodcounts1$neutropenia, useNA = "ifany") # 15861
summary(covs_bloodcounts1$Neutrophill.count)  # 15861 NAs

table(covs_bloodcounts1$cytopenia, useNA = "ifany")


## get tobacco use disorder
tobacc = fread(
    "/medpop/esp2/mzekavat/UKBB/PhenoFiles/PheCODEs_new/2021-01-08_ukb_phecode_MentalNeuroSensorySymptoms_March2020fu.csv", 
    data.table = F)
tobacc_use_disorder <- tobacc %>%
    select(id = f.eid, Tobacco_use_disorder_ANY)


covs_final <- left_join(covs, tobacc_use_disorder, by="id") %>%
    left_join(covs_bloodcounts1 %>% select(id, cytopenia), by ="id")


write.table(covs_final, "Data/ukb_phewas/ukb_phewas_covs_final.tsv", sep="\t", row.names = F, quote = F)
