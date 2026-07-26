###### 01 - library #####
library(dplyr)
library(data.table)

gc()
rm(list=ls())

setwd("/medpop/esp2/lli/chip_gwas_rev/")


###### 02 - files #####
# covariates -------------------------------------------------------------------
covs = fread("/medpop/esp2/mzekavat/Partners_BB/phenos/all_basic_MGBB_phenos.plus_mCAs_CHIP.txt", data.table = F)


## note that the filtering command from dplyr and which command from Maryam's old code is not the same - 
## dplyr generates a specific NA claude while 'which' treats all the NAs as FALSE
covs_clean <- covs %>%
    filter(!(!is.na(sexCheck_toRemove) & sexCheck_toRemove == 1)) %>%
    filter(!((!is.na(Prev_Heme_CA) & Prev_Heme_CA == 1)|(!is.na(age)&age< 20))) %>%
    mutate(ever_smoked = factor(Smoking_Final)) %>%
    select(id = Biobank.Subject.ID, Prev_Heme_CA, age, age2, Race, ever_smoked, sex=Gender)

table(covs$Smoking_Final, useNA = "always") # no NA
table(covs$Race, useNA = "always") # no NA


knn_df <- fread("/medpop/esp2/projects/MGB_Biobank/genotype/53K_GSA/release/GSA_53K.hg38.1KGPCA.knn_rf.tsv.gz", data.table = F)
knn_df <- knn_df %>%
    select(IID, paste("PC",1:10, sep=""), knn) %>%
    mutate(id = sub(".*-", "", IID), .after=IID) %>%
    select(id, paste("PC",1:10,sep=""), knn) %>%
    mutate(id = as.numeric(id)) %>%
    filter(!is.na(id))


merge <- left_join(covs_clean, knn_df, by = "id")
table(merge$Race, merge$knn, useNA = "ifany")

tobacc = fread(
    "/medpop/esp2/mzekavat/Partners_BB/phenos/PheCodes/clean_mgbb.ICDdiagnosis.phecodes2020.MentalNeuroSensorySymptoms.20201105_gender.restricted.csv", 
    data.table=F)

tobacc_df <- tobacc %>%
    select(id = Biobank.Subject.ID, Tobacco_use_disorder_ANY)
rm(tobacc)

# I think the best way, since we are stratifying by knn, is to control for PC1-10 and drop Race
covs_final <- merge %>%
    left_join(tobacc_df, by = "id") %>%
    select(id, Prev_Heme_CA, age, age2, ever_smoked, contains("PC"), knn, sex, Tobacco_use_disorder_ANY)

write.table(covs_final, "Data/mgb_phewas/mgbb_covs_final.tsv", sep='\t', row.names = F, quote = F)




# bloodcounts -------------------------------------------------------------
gc()
rm(list = ls())

ab_raw <- fread("Data/mgb_phewas/bloodcounts/data_all_blood.txt", data.table = F)
ab_raw <- ab_raw %>%
    rename(id = Subject_Id)
covs_final <- fread("Data/mgb_phewas/mgbb_covs_final.tsv", data.table = F)


ab_covs <- left_join(ab_raw, covs_final %>% select(id, sex), by = "id")
ab_covs1 <- ab_covs %>%
    mutate(
        anemia = case_when(
            is.na(HGB) ~ NA,
            (sex == "F" & HGB < 12) | (sex == "M" & HGB < 13) ~ 1,
            (sex == "F" & HGB >= 12) | (sex == "M" & HGB >= 13) ~ 0,
            HGB >= 13 ~ 0,
            HGB < 12  ~ 1,
            TRUE ~ NA
        ),
        thrombocytopenia = case_when(
            is.na(PLT) ~ NA,
            PLT < 150 ~ 1,
            PLT >= 150 ~ 0,
            TRUE ~ NA
        ),
        any_cytopenia = case_when(
            is.na(anemia) & is.na(thrombocytopenia) ~ NA,
            anemia == 1 | thrombocytopenia == 1 ~ 1,
            TRUE ~ 0
        )
    )


table(ab_covs1$any_cytopenia, useNA = "ifany")
# 0       1 
# 21296  7763


write.table(ab_covs1 %>% select(id, any_cytopenia), "Data/mgb_phewas/mgbb_cytopenia_final.tsv",
            sep="\t", row.names = F, quote = F)





### all of the ukbb phecodes data
all_phecodes <- list.files("/medpop/esp2/mzekavat/Partners_BB/phenos/PheCodes/")
phecodes_cat <- all_phecodes[!grepl("OLD", all_phecodes)]
# "/medpop/esp2/mzekavat/UKBB/PhenoFiles/PheCODEs_new/2021-01-08_ukb_phecode_incident_March2020fu.csv" 
# is coded binary for incidence cases
# incident.days is for number of days till incidence
# 2021-01-08_ukb_phecode_any_March2020fu.csv also binary

phecodes_cat_final <- gsub("clean_mgbb.ICDdiagnosis.phecodes2020.|.20201105_gender.restricted.csv", "", phecodes_cat)

write.table(phecodes_cat_final, "Data/mgb_phewas/phecode_list.tsv", sep="\t", row.names = F, quote = F, col.names = F)











